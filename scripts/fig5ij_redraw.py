"""Fig5I / Fig5J — Wet experiment validation bar plots.

目标：与 Fig5H 重做版风格严格统一.
  - 配色 PAL = {"Control":#7F7F7F, "PRRSV":#D62728, "APRR":#1F77B4}
  - Arial 字体; svg.fonttype='none' 使文字保留为 <text> 元素 (AI 可编辑)
  - 去掉上 / 右脊柱, 左 / 下深灰 #2A2A2A
  - 误差线 SD (与用户原图同宽感, 来自汇总表 sd 列)
  - *** 仅画 Control vs PRRSV 与 PRRSV vs APRR (按 Welch 近似 t 检验)
  - "down X.X%" 蓝色加粗文字置于 APRR 柱顶 (X.X 由 (PRRSV-APRR)/PRRSV × 100)

数据源:
  - 肺 p-p65 IF: wet_experiment_archive/data/02_summary_pp65_lung_IF.csv (n=4/组)
  - 血清 IL-1β: wet_experiment_archive/data/01_summary_IL1beta_serum_ELISA.csv (n=6/组)
约束:
  - 唯一治疗组 = APRR (禁用剂量标签). CSV 中 "APRR4" 在 X tick 上映射为 "APRR"
  - 所有可交付 SVG 保留 <text> 元素 (不转曲)
"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

warnings.filterwarnings("ignore")

DATA_DIR = (
    r"data/"
    r"scRNA-seq/reanalysis_v2/wet_experiment_archive/data"
)
FIG5_ROOT = (
    r"data/"
    r"scRNA-seq/reanalysis_v2/figures_final/Fig5"
)
DIR_VAL = os.path.join(FIG5_ROOT, "validation")
DIR_MSS = os.path.join(FIG5_ROOT, "Fig5_Manuscript")

PAL = {"Control": "#7F7F7F", "PRRSV": "#D62728", "APRR": "#1F77B4"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
})


def welch_ttest_from_summary(m1, s1, n1, m2, s2, n2):
    """Welch's two-sample t from summary stats → p (two-sided)."""
    se1 = s1 ** 2 / n1
    se2 = s2 ** 2 / n2
    se = np.sqrt(se1 + se2)
    if se == 0:
        return 1.0
    t = (m1 - m2) / se
    df = (se1 + se2) ** 2 / (se1 ** 2 / (n1 - 1) + se2 ** 2 / (n2 - 1))
    return 2 * stats.t.sf(np.abs(t), df)


def style_ax(ax, ymax, ylabel, y_locator):
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["Control", "PRRSV", "APRR"], rotation=30, ha="right")
    ax.set_xlim(0.4, 3.6)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(ylabel, labelpad=4)
    ax.yaxis.set_major_locator(plt.MultipleLocator(y_locator))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#2A2A2A")
    ax.tick_params(axis="both", colors="#2A2A2A")


def draw_bar(ax, m, sd, position, color, width=0.74, cap=0.10):
    """Filled bar + vertical SD error bar with end caps (no edge)."""
    ax.bar(position, m, width=width, color=color, edgecolor="none",
           alpha=0.95, zorder=2)
    # error bar
    ax.plot([position, position], [m - sd, m + sd],
            color="#1A1A1A", lw=0.8, zorder=3, solid_capstyle="butt")
    ax.plot([position - cap, position + cap], [m + sd, m + sd],
            color="#1A1A1A", lw=0.8, zorder=3, solid_capstyle="butt")
    ax.plot([position - cap, position + cap], [m - sd, m - sd],
            color="#1A1A1A", lw=0.8, zorder=3, solid_capstyle="butt")


def ann_p(ax, x1, x2, y_top, h, p, group_order):
    """Draw *** bar above two columns if p<0.001 (this paper's standard)."""
    if p >= 0.001:
        return
    ax.plot([x1, x1, x2, x2], [y_top, y_top + h, y_top + h, y_top],
            color="#222222", lw=0.6)
    ax.text((x1 + x2) / 2.0, y_top + h + h * 0.10, "***",
            ha="center", va="bottom", fontsize=8.5, color="#222222",
            fontweight="bold")


def add_downlabel(ax, m_prrsv, m_aprr, position_aprr, color="#1F77B4"):
    """'down X.X%' label above the APRR bar (blue, bold)."""
    pct = (m_prrsv - m_aprr) / m_prrsv * 100
    xmid = position_aprr
    ytop = m_aprr + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.020
    # 0.04 of y-range above top
    ypos = ax.get_ylim()[1] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.06
    ax.text(xmid, ypos, f"down {pct:.1f}%",
            ha="center", va="top",
            fontsize=9, color=color, fontweight="bold")


def render_panel(df, name_suffix, title, ylabel, y_locator,
                 panel_id, down_label=True):
    """Single-panel bar plot — used for both Fig5I (pp65) and Fig5J (IL-1β)."""
    fig = plt.figure(figsize=(3.4, 3.4), dpi=300)
    gs = fig.add_gridspec(1, 1, left=0.20, right=0.96, top=0.86, bottom=0.28)
    ax = fig.add_subplot(gs[0, 0])

    # df already has rows 'Control','PRRSV','APRR' in order
    rows = df.set_index("group").loc[["Control", "PRRSV", "APRR"]]
    means = rows["mean"].values
    sds = rows["sd"].values
    ns = rows["n"].values
    colors = [PAL["Control"], PAL["PRRSV"], PAL["APRR"]]

    # ymax — 给 PRRSV 顶 + SD (含 *** 高度) + down X% 标签充足顶部空隙
    # reserve = max(SD_PRRSV) + 0.16 * (PRRSV+SD) height for two *** bars + "down" label
    span = means.max() + sds.max()
    ymax = float(np.ceil(span * 1.45 / 10) * 10)

    style_ax(ax, ymax, ylabel, y_locator)
    ax.set_title(title, pad=5)

    for x, m, sd, c in zip([1, 2, 3], means, sds, colors):
        draw_bar(ax, m, sd, x, c)

    # p values: Ctrl vs PRRSV, PRRSV vs APRR (Welch summary t)
    m_c, s_c, n_c = rows.loc["Control", ["mean", "sd", "n"]]
    m_p, s_p, n_p = rows.loc["PRRSV", ["mean", "sd", "n"]]
    m_a, s_a, n_a = rows.loc["APRR", ["mean", "sd", "n"]]

    p_cp = welch_ttest_from_summary(m_c, s_c, n_c, m_p, s_p, n_p)
    p_pa = welch_ttest_from_summary(m_p, s_p, n_p, m_a, s_a, n_a)

    # *** bars (上层 PRRSV-vs-APRR, 下层 Ctrl-vs-PRRSV, 互不重叠)
    h = ymax * 0.025                  # *** 横线短腿高
    bar_y_low = ymax * 0.80           # Ctrl vs PRRSV
    bar_y_hi  = ymax * 0.86           # PRRSV vs APRR
    ann_p(ax, 1, 2, bar_y_low, h, p_cp, ["Control", "PRRSV"])
    ann_p(ax, 2, 3, bar_y_hi,  h, p_pa, ["PRRSV", "APRR"])

    if down_label:
        # 'down X.X%' 蓝色加粗 — 放在 APRR 柱内顶上方单独一行 (错开 ***)
        pct = (m_p - m_a) / m_p * 100
        ypos = ymax * 0.99            # 顶部留 1% 余量
        ax.text(3, ypos, f"down {pct:.1f}%",
                ha="center", va="top",
                fontsize=8.5, color="#1F77B4", fontweight="bold")

    out_base = f"Fig5{panel_id}_"
    fname = out_base + name_suffix
    out_png = os.path.join(DIR_VAL, fname + ".png")
    out_svg = os.path.join(DIR_VAL, fname + ".svg")
    out_pdf = os.path.join(DIR_VAL, fname + ".pdf")
    fig.savefig(out_png, dpi=600, facecolor="white")
    fig.savefig(out_svg, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)

    return {
        "name": fname,
        "png": out_png,
        "svg": out_svg,
        "pdf": out_pdf,
        "p_cp": p_cp,
        "p_pa": p_pa,
        "down_pct": (m_p - m_a) / m_p * 100,
        "means": (m_c, m_p, m_a),
        "sds": (s_c, s_p, s_a),
        "ns": (int(n_c), int(n_p), int(n_a)),
    }


def mirror_to_manuscript(fnames):
    """Copy PNG/SVG/PDF to Fig5_Manuscript (canonical)."""
    import shutil
    for f in fnames:
        for ext in (".png", ".svg", ".pdf"):
            src = os.path.join(DIR_VAL, f + ext)
            if not os.path.exists(src):
                continue
            dst = os.path.join(DIR_MSS, f + ext)
            shutil.copyfile(src, dst)


# ============== Load data (group CSV label "APRR4" → "APRR") ==============
df_i = pd.read_csv(os.path.join(DATA_DIR, "02_summary_pp65_lung_IF.csv"))
df_i["group"] = df_i["group"].replace({"APRR4": "APRR"})

df_j = pd.read_csv(os.path.join(DATA_DIR, "01_summary_IL1beta_serum_ELISA.csv"))
df_j["group"] = df_j["group"].replace({"APRR4": "APRR"})

# ============== Fig5I ==============
out_i = render_panel(
    df_i,
    name_suffix="pp65_IF_bar",
    title="Lung p-p65 (NF-κB activation)",
    ylabel="Relative intensity (a.u.)",
    y_locator=20,
    panel_id="I",
)

# ============== Fig5J ==============
out_j = render_panel(
    df_j,
    name_suffix="IL1b_ELISA_bar",
    title="Serum IL-1β (ELISA)",
    ylabel="IL-1β (pg/mL)",
    y_locator=50,
    panel_id="J",
)

# ============== Mirror ==============
mirror_to_manuscript([out_i["name"], out_j["name"]])

# ============== Stats summary ==============
print("==== Fig5I (p-p65 IF) ====")
print(f"  n    = {out_i['ns']}")
print(f"  mean = {[f'{m:.2f}' for m in out_i['means']]}")
print(f"  sd   = {[f'{s:.2f}' for s in out_i['sds']]}")
print(f"  Ctrl vs PRRSV p = {out_i['p_cp']:.3e}")
print(f"  PRRSV vs APRR p = {out_i['p_pa']:.3e}")
print(f"  PRRSV→APRR down = {out_i['down_pct']:.1f}%")
print()
print("==== Fig5J (IL-1β ELISA) ====")
print(f"  n    = {out_j['ns']}")
print(f"  mean = {[f'{m:.2f}' for m in out_j['means']]}")
print(f"  sd   = {[f'{s:.2f}' for s in out_j['sds']]}")
print(f"  Ctrl vs PRRSV p = {out_j['p_cp']:.3e}")
print(f"  PRRSV vs APRR p = {out_j['p_pa']:.3e}")
print(f"  PRRSV→APRR down = {out_j['down_pct']:.1f}%")
