# -*- coding: utf-8 -*-
"""
build_wet_archive.py
====================
将 APRR 抗 PRRSV 项目的全部湿实验结果统一排版、统计、可视化，并归档到
reanalysis_v2/wet_experiment_archive/ 下。

产物：
  wet_experiment_archive/
    ├── data/                     排版后的规范汇总表（含 mean/SD/SEM/n）
    │   ├── 01_summary_IL1beta_serum_ELISA.csv
    │   ├── 02_summary_pp65_lung_IF.csv
    │   ├── 03_summary_pp65_image_quant.csv
    │   ├── 04_summary_viral_load_by_tissue.csv
    │   ├── 05_summary_viral_load_serum_timeseries.csv
    │   ├── 06_summary_phenotype.csv
    │   ├── 07_summary_antiviral_dose.csv
    │   └── raw/                  原始 wet CSV 原样复制
    ├── figures/                  统计图表（PNG，300 dpi）
    │   ├── fig_IL1beta_serum_ELISA.png
    │   ├── fig_pp65_lung_IF.png
    │   ├── fig_pp65_image_quant.png
    │   ├── fig_viral_load_serum.png
    │   ├── fig_viral_load_by_tissue.png
    │   ├── fig_phenotype_temperature.png
    │   ├── fig_phenotype_weight.png
    │   ├── fig_phenotype_clinical_score.png
    │   ├── fig_phenotype_survival.png
    │   └── fig_antiviral_dose.png
    └── README.md                 归档索引（由单独步骤生成）

分组约定：
  - 大部分湿实验用 APRR4（与 PRRSV 同期给药组）
  - p-p65 肺 IF 用户原始组名为 APRR，与 APRR4 为同批动物不同检测命名；
    本脚本在汇总表中统一标记为 APRR(=APRR4) 以便阅读，并在 README 标注。

统计方法：
  - 连续变量以 mean ± SD 报告，n 为动物数/重复数。
  - 组间比较：PRRSV vs APRR(4) 独立样本 t 检验（双尾），
    显著性 *p<0.05, **p<0.01, ***p<0.001。
"""

import os
import shutil
import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 注册系统中文黑体（SimHei），避免图中中文标签乱码
_CN_FONT_CANDIDATES = [
    r"C:\Windows\Fonts\simhei.ttf",
    r"C:\Windows\Fonts\msyh.ttc",
    r"C:\Windows\Fonts\simsun.ttc",
]
_CN_FONT = None
for _fp in _CN_FONT_CANDIDATES:
    if os.path.exists(_fp):
        try:
            font_manager.fontManager.addfont(_fp)
            _CN_FONT = font_manager.FontProperties(fname=_fp).get_name()
            break
        except Exception:
            continue
if _CN_FONT:
    plt.rcParams["font.family"] = _CN_FONT
    plt.rcParams["axes.unicode_minus"] = False
    print(f"[OK] 中文字体已注册: {_CN_FONT}")
else:
    print("[WARN] 未找到系统中文字体，中文标签可能显示为方框")

# ----------------------------------------------------------------------------
# 路径配置
# ----------------------------------------------------------------------------
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # reanalysis_v2
WET = os.path.join(BASE, "tables", "wet")
OUT = os.path.join(BASE, "wet_experiment_archive")
DATA = os.path.join(OUT, "data")
FIGS = os.path.join(OUT, "figures")
RAW = os.path.join(DATA, "raw")

os.makedirs(DATA, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)
os.makedirs(RAW, exist_ok=True)

GROUP_ORDER = ["Control", "PRRSV", "APRR4"]

# 组织名中英映射（用于图表 axis，避免纯中文依赖字体）
TISSUE_EN = {
    "心脏": "Heart", "肝脏": "Liver", "脾脏": "Spleen", "肺脏": "Lung",
    "肾脏": "Kidney", "空肠": "Jejunum", "结肠": "Colon", "脑": "Brain",
    "扁桃体Tonsil": "Tonsil", "下颌淋巴Submaxillary nodes": "Submaxillary LN",
    "腹股沟淋巴Inguinal nodes": "Inguinal LN",
}
# 指标名中英映射
METRIC_EN = {
    "体温变化改改": "Body temperature", "体重变化改改": "Body weight change",
    "临床症状评分": "Clinical score", "生存曲线": "Survival",
}

def en_tissue(t):
    return TISSUE_EN.get(t, t)

def en_metric(m):
    return METRIC_EN.get(m, m)

# 统一绘图样式
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

COLORS = {"Control": "#4C72B0", "PRRSV": "#C44E52", "APRR4": "#55A868"}


def sig_label(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def summarize_long(df, group_col, value_col):
    """对长格式数据按 group 计算 mean/sd/sem/n。"""
    rows = []
    for g in GROUP_ORDER:
        sub = df[df[group_col] == g][value_col].dropna().values.astype(float)
        if len(sub) == 0:
            continue
        rows.append({
            "group": g,
            "n": len(sub),
            "mean": round(float(np.mean(sub)), 3),
            "sd": round(float(np.std(sub, ddof=1)), 3) if len(sub) > 1 else 0.0,
            "sem": round(float(np.std(sub, ddof=1) / np.sqrt(len(sub))), 3) if len(sub) > 1 else 0.0,
        })
    return pd.DataFrame(rows)


def ttest_prrsv_vs_aprr(df, group_col, value_col):
    """返回 (t, p) 比较 PRRSV vs APRR4(APRR)。"""
    a = df[df[group_col] == "PRRSV"][value_col].dropna().values.astype(float)
    b = df[df[group_col].isin(["APRR4", "APRR"])][value_col].dropna().values.astype(float)
    if len(a) < 2 or len(b) < 2:
        return (np.nan, np.nan)
    return stats.ttest_ind(a, b, equal_var=False)


def bar_with_sig(ax, summary, title, ylabel, compare=True, fmt="{:.1f}"):
    groups = summary["group"].tolist()
    means = summary["mean"].values
    sds = summary["sd"].values
    xs = np.arange(len(groups))
    bars = ax.bar(xs, means, yerr=sds, capsize=6,
                  color=[COLORS[g] for g in groups], edgecolor="black", alpha=0.9)
    ax.set_xticks(xs)
    ax.set_xticklabels(groups)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    for x, m in zip(xs, means):
        ax.text(x, m, fmt.format(m), ha="center", va="bottom", fontsize=10)
    if compare and "PRRSV" in groups and ("APRR4" in groups or "APRR" in groups):
        i_p = groups.index("PRRSV")
        i_a = groups.index("APRR4") if "APRR4" in groups else groups.index("APRR")
        ymax = max(means[i_p] + sds[i_p], means[i_a] + sds[i_a])
        ytop = ymax * 1.12
        ax.plot([i_p, i_p, i_a, i_a], [ymax * 1.02, ytop, ytop, ymax * 1.02],
                lw=1.2, color="black")
        # p 值由调用方通过 summary.attrs 传入
        p = summary.attrs.get("p_prrsv_vs_aprr", np.nan) if hasattr(summary, "attrs") else np.nan
        lab = sig_label(p) if not (isinstance(p, float) and np.isnan(p)) else ""
        ax.text((i_p + i_a) / 2, ytop, lab, ha="center", va="bottom", fontsize=12)
    return ax


# ============================================================================
# 1. IL-1beta 血清 ELISA
# ============================================================================
def build_il1beta():
    df = pd.read_csv(os.path.join(WET, "wet_serum_cytokines.csv"))
    df = df.rename(columns={"value_pg_mL": "value"})
    summary = summarize_long(df, "group", "value")
    summary.attrs["p_prrsv_vs_aprr"] = ttest_prrsv_vs_aprr(df, "group", "value")[1]
    # 倍数变化
    m = summary.set_index("group")["mean"]
    fc_virus = m["PRRSV"] / m["Control"]
    fc_treat = m["APRR4"] / m["PRRSV"]
    summary["fold_change_vs_Control"] = [
        "-" if g == "Control" else round(m[g] / m["Control"], 3) for g in summary["group"]
    ]
    summary["fold_change_vs_PRRSV"] = [
        "-" if g in ("Control", "PRRSV") else round(m[g] / m["PRRSV"], 3) for g in summary["group"]
    ]
    out = os.path.join(DATA, "01_summary_IL1beta_serum_ELISA.csv")
    summary.to_csv(out, index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bar_with_sig(ax, summary, "Serum IL-1β (ELISA)", "IL-1β (pg/mL)")
    fig.savefig(os.path.join(FIGS, "fig_IL1beta_serum_ELISA.png"))
    plt.close(fig)
    print(f"[OK] IL-1β ELISA -> {out}  (PRRSV/Control={fc_virus:.2f}x, APRR4/PRRSV={fc_treat:.2f}x)")
    return summary


# ============================================================================
# 2. p-p65 肺 IF（用户实测）
# ============================================================================
def build_pp65_if():
    df = pd.read_csv(os.path.join(WET, "wet_lung_pp65_if.csv"))
    df = df.rename(columns={"value": "value"})
    # 组名 APRR -> APRR4 以便统一阅读
    df["group"] = df["group"].replace({"APRR": "APRR4"})
    summary = summarize_long(df, "group", "value")
    summary.attrs["p_prrsv_vs_aprr"] = ttest_prrsv_vs_aprr(df, "group", "value")[1]
    m = summary.set_index("group")["mean"]
    fc_virus = m["PRRSV"] / m["Control"]
    fc_treat = m["APRR4"] / m["PRRSV"]
    summary["fold_change_vs_Control"] = [
        "-" if g == "Control" else round(m[g] / m["Control"], 3) for g in summary["group"]
    ]
    summary["fold_change_vs_PRRSV"] = [
        "-" if g in ("Control", "PRRSV") else round(m[g] / m["PRRSV"], 3) for g in summary["group"]
    ]
    out = os.path.join(DATA, "02_summary_pp65_lung_IF.csv")
    summary.to_csv(out, index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bar_with_sig(ax, summary, "Lung p-p65 (Immunofluorescence)", "p-p65 IF index")
    fig.savefig(os.path.join(FIGS, "fig_pp65_lung_IF.png"))
    plt.close(fig)
    print(f"[OK] p-p65 IF -> {out}  (PRRSV/Control={fc_virus:.2f}x, APRR4/PRRSV={fc_treat:.2f}x)")
    return summary


# ============================================================================
# 3. p-p65 图像定量（脚本产出，4 子区域重复）
# ============================================================================
def build_pp65_image_quant():
    qdir = os.path.join(WET, "quant_pp65")
    src = os.path.join(qdir, "pp65_grid_summary.csv")
    df = pd.read_csv(src)
    # 该 summary 已有 rep1-4, mean, sd, n_valid
    out = os.path.join(DATA, "03_summary_pp65_image_quant.csv")
    df.to_csv(out, index=False, encoding="utf-8-sig")
    # 复制 detail 一并保留
    shutil.copy(os.path.join(qdir, "pp65_grid_detail.csv"),
                os.path.join(DATA, "03_detail_pp65_image_quant.csv"))

    # 点图：每张图 4 子区域重复
    detail = pd.read_csv(os.path.join(qdir, "pp65_grid_detail.csv"))
    if {"file", "value"}.issubset(detail.columns):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        files = detail["file"].unique()
        xs = np.arange(len(files))
        for i, f in enumerate(files):
            vals = detail[detail["file"] == f]["value"].values.astype(float)
            ax.scatter(np.full_like(vals, i) + np.random.uniform(-0.12, 0.12, len(vals)),
                       vals, color="#333333", alpha=0.7, zorder=3)
            ax.scatter([i], [np.mean(vals)], color="#C44E52", s=60, marker="D", zorder=4)
        ax.set_xticks(xs)
        ax.set_xticklabels([os.path.splitext(f)[0] for f in files], rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("Mean R-channel intensity (4 sub-regions)")
        ax.set_title("p-p65 image quantitation (4 sub-region replicates)")
        ax.set_xlabel("Image (red single-channel, 4 sub-region reps)")
        fig.savefig(os.path.join(FIGS, "fig_pp65_image_quant.png"))
        plt.close(fig)
        print(f"[OK] p-p65 image quant -> {out}  ({len(files)} images x 4 sub-regions)")
    else:
        print(f"[WARN] p-p65 detail 列名不符，跳过点图: {detail.columns.tolist()}")
    return df


# ============================================================================
# 4. 病毒载量（按组织）+ 5. 血清时间序列
# ============================================================================
def build_viral_load():
    df = pd.read_csv(os.path.join(WET, "wet_viral_load_by_tissue.csv"))

    # 4) 各组织汇总（排除聚合项 血清/器官/3day/6day/10day 血清，单独处理）
    agg_terms = {"血清", "器官", "3day 血清", "6day 血清", "10day 血清"}
    tissues = [t for t in df["tissue"].unique() if t not in agg_terms]
    rows = []
    for t in tissues:
        sub = df[df["tissue"] == t]
        s = summarize_long(sub, "group", "value")
        s.insert(0, "tissue", t)
        rows.append(s)
    vload = pd.concat(rows, ignore_index=True)
    out = os.path.join(DATA, "04_summary_viral_load_by_tissue.csv")
    vload.to_csv(out, index=False, encoding="utf-8-sig")

    # 5) 血清时间序列（3day/6day/10day + 总体血清）
    series_terms = {"3day 血清": "3 dpc", "6day 血清": "6 dpc", "10day 血清": "10 dpc", "血清": "All"}
    srows = []
    for t, label in series_terms.items():
        sub = df[df["tissue"] == t]
        s = summarize_long(sub, "group", "value")
        s.insert(0, "timepoint", label)
        srows.append(s)
    vload_series = pd.concat(srows, ignore_index=True)
    out2 = os.path.join(DATA, "05_summary_viral_load_serum_timeseries.csv")
    vload_series.to_csv(out2, index=False, encoding="utf-8-sig")

    # 图：总体血清病毒载量
    serum = df[df["tissue"] == "血清"]
    serum_s = summarize_long(serum, "group", "value")
    serum_s.attrs["p_prrsv_vs_aprr"] = ttest_prrsv_vs_aprr(serum, "group", "value")[1]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bar_with_sig(ax, serum_s, "Serum viral load", "Log10 viral copies (a.u.)")
    fig.savefig(os.path.join(FIGS, "fig_viral_load_serum.png"))
    plt.close(fig)

    # 图：各组织病毒载量（Control 全 0，仅画 PRRSV vs APRR4）
    fig, ax = plt.subplots(figsize=(9, 5))
    means = []
    sds = []
    labels = []
    for t in tissues:
        sub = df[df["tissue"] == t]
        pr = sub[sub["group"] == "PRRSV"]["value"].values.astype(float)
        ap = sub[sub["group"] == "APRR4"]["value"].values.astype(float)
        means.append([np.mean(pr), np.mean(ap)])
        sds.append([np.std(pr, ddof=1), np.std(ap, ddof=1)])
        labels.append(t)
    means = np.array(means)
    sds = np.array(sds)
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w/2, means[:, 0], yerr=sds[:, 0], width=w, capsize=4,
           color=COLORS["PRRSV"], edgecolor="black", label="PRRSV")
    ax.bar(x + w/2, means[:, 1], yerr=sds[:, 1], width=w, capsize=4,
           color=COLORS["APRR4"], edgecolor="black", label="APRR4")
    ax.set_xticks(x)
    ax.set_xticklabels([en_tissue(t) for t in labels], rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Log10 viral copies (a.u.)")
    ax.set_title("Tissue viral load: PRRSV vs APRR4")
    ax.legend()
    fig.savefig(os.path.join(FIGS, "fig_viral_load_by_tissue.png"))
    plt.close(fig)
    print(f"[OK] viral load -> {out}, {out2}")
    return vload, vload_series


# ============================================================================
# 6. 表型（体温/体重/临床/生存）
# ============================================================================
def build_phenotype():
    df = pd.read_csv(os.path.join(WET, "wet_phenotype.csv"))
    metrics = df["metric"].unique()
    rows = []
    for mt in metrics:
        sub = df[df["metric"] == mt]
        s = summarize_long(sub, "group", "value")
        s.insert(0, "metric", mt)
        rows.append(s)
    pheno = pd.concat(rows, ignore_index=True)
    out = os.path.join(DATA, "06_summary_phenotype.csv")
    pheno.to_csv(out, index=False, encoding="utf-8-sig")

    # 体温 / 体重 / 临床评分：取末次观测或用 mean 作图（此处按组 mean±sd 概览）
    metric_map = {
        "体温变化改改": ("Body temperature (°C)", "fig_phenotype_temperature.png"),
        "体重变化改改": ("Body weight change (%)", "fig_phenotype_weight.png"),
        "临床症状评分": ("Clinical score", "fig_phenotype_clinical_score.png"),
    }
    for mt, (ylabel, fname) in metric_map.items():
        sub = df[df["metric"] == mt]
        s = summarize_long(sub, "group", "value")
        if s.empty:
            continue
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        bar_with_sig(ax, s, en_metric(mt), ylabel, compare=False)
        fig.savefig(os.path.join(FIGS, fname))
        plt.close(fig)

    # 生存曲线：Kaplan-Meier
    surv = df[df["metric"] == "生存曲线"]
    if not surv.empty:
        from scipy.interpolate import interp1d
        fig, ax = plt.subplots(figsize=(6, 4.5))
        for g in GROUP_ORDER:
            sub = surv[surv["group"] == g]
            if sub.empty:
                continue
            vals = np.sort(sub["value"].values.astype(float))
            # value 假定为生存天数；绘制经验生存函数
            n = len(vals)
            surv_frac = np.arange(1, n + 1)[::-1] / n
            ax.step(np.concatenate([[0], vals]), np.concatenate([[1.0], surv_frac]),
                    where="post", label=g, color=COLORS[g], lw=2)
        ax.set_xlabel("Days")
        ax.set_ylabel("Survival fraction")
        ax.set_title("Survival curve")
        ax.set_ylim(0, 1.05)
        ax.legend()
        fig.savefig(os.path.join(FIGS, "fig_phenotype_survival.png"))
        plt.close(fig)
        print(f"[OK] phenotype -> {out} (含生存曲线)")
    else:
        print(f"[OK] phenotype -> {out}")
    return pheno


# ============================================================================
# 7. 抗病毒剂量效应（细胞实验）
# ============================================================================
def build_dose():
    df = pd.read_csv(os.path.join(WET, "wet_antiviral_dose.csv"))
    rows = []
    for a in df["assay"].unique():
        sub = df[df["assay"] == a]
        s = sub.groupby("dose")["value"].agg(["mean", "std", "count"]).reset_index()
        s = s.rename(columns={"mean": "mean", "std": "sd", "count": "n"})
        s.insert(0, "assay", a)
        s["dose"] = s["dose"].astype(int)
        rows.append(s)
    dose = pd.concat(rows, ignore_index=True)
    out = os.path.join(DATA, "07_summary_antiviral_dose.csv")
    dose.to_csv(out, index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for a in df["assay"].unique():
        sub = df[df["assay"] == a]
        s = sub.groupby("dose")["value"].mean()
        ax.plot(s.index.astype(int), s.values, marker="o", label=a)
    ax.set_xlabel("Dose (a.u.)")
    ax.set_ylabel("Antiviral effect")
    ax.set_title("APRR dose-response (cell assay)")
    ax.legend(title="Assay")
    fig.savefig(os.path.join(FIGS, "fig_antiviral_dose.png"))
    plt.close(fig)
    print(f"[OK] antiviral dose -> {out}")
    return dose


# ============================================================================
# 复制原始 wet 文件
# ============================================================================
def copy_raw():
    copied = []
    for f in os.listdir(WET):
        src = os.path.join(WET, f)
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(RAW, f))
            copied.append(f)
    # quant_pp65 子目录整体复制
    qsrc = os.path.join(WET, "quant_pp65")
    if os.path.isdir(qsrc):
        qdst = os.path.join(RAW, "quant_pp65")
        os.makedirs(qdst, exist_ok=True)
        for f in os.listdir(qsrc):
            shutil.copy(os.path.join(qsrc, f), os.path.join(qdst, f))
        copied.append("quant_pp65/")
    print(f"[OK] raw copied ({len(copied)} items): {copied}")


def main():
    print("=== Building wet experiment archive ===")
    build_il1beta()
    build_pp65_if()
    build_pp65_image_quant()
    build_viral_load()
    build_phenotype()
    build_dose()
    copy_raw()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
