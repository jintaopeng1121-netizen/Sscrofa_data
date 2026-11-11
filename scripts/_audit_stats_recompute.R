## 复算稿件体内所有统计量，判定 p 值与检验方法是否一致
## 输出：tables/_audit_stats_recompute.txt（人工核对用，不进投稿）
suppressMessages({ library(Seurat) })
OUT <- "data/scRNA-seq/reanalysis_v2"

rep_row <- function(label, ctrl, prrsv, aprr, n_ctrl=NULL) {
  a <- tryCatch(t.test(aprr, prrsv, var.equal = TRUE),  error=function(e) NULL)
  b <- tryCatch(t.test(aprr, prrsv, var.equal = FALSE), error=function(e) NULL)
  cat(sprintf("%-22s  mean APRR=%.3f  PRRSV=%.3f  diff=%+.3f\n", label,
              mean(aprr), mean(prrsv), mean(aprr)-mean(prrsv)))
  cat(sprintf("   Student (equal var): t=%.3f df=%.2f p=%.4g\n", a$statistic, a$parameter, a$p.value))
  cat(sprintf("   Welch              : t=%.3f df=%.2f p=%.4g\n", b$statistic, b$parameter, b$p.value))
  cat(sprintf("   SD: APRR=%.3f PRRSV=%.3f | n: APRR=%d PRRSV=%d\n\n",
              sd(aprr), sd(prrsv), length(aprr), length(prrsv)))
}

cat("=== 1. 组织/血清病毒载量（log10 copies）===\n")
rep_row("lung",  c(0,0,0),            c(9.4,9.09,9.12),        c(7.56,7.7,8.32))
rep_row("spleen",c(0,0,0),            c(8.05,8.69,7.78),      c(7.2,7.57,7.29))
rep_row("serum 10dpi", c(0,0,0,0),    c(8.43,8.89,7.03,7.39), c(6.48,7.19,7.19,7.47))

cat("=== 2. 血清 IL-1beta (pg/mL) ===\n")
rep_row("IL-1beta", c(56.74,49.72,55.92,65.38,56.27,87.17),
        c(179.39,187.66,149.58,129.31,202.49,212.17),
        c(79.52,69.23,89.62,86.98,96.49,113.07))
cat("   reduction = %.1f%%\n\n", (1-mean(c(79.52,69.23,89.62,86.98,96.49,113.07))/mean(c(179.39,187.66,149.58,129.31,202.49,212.17)))*100)

cat("=== 3. 临床评分 / 体温 / 体重（原始为重复测量扁平数组）===\n")
ctrl_sc <- c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
prr_sc  <- c(0,1,4,3,5,6,9,9,10,10,0,0,1,4,6,8,6,7,8,7,9,0,4,4,2,5,6,7,6,8,7,8,0,2,3,4,2,7,7,5,7,8,6,0,0,3,4,5,3,4,4,6,8,6,0,2,0,6,7,7,6,7,5)
apr_sc  <- c(0,1,4,3,2,7,7,6,7,5,4,0,0,1,1,1,3,6,7,4,4,5,0,2,1,1,5,4,5,5,5,7,5,0,0,1,3,3,6,4,3,4,3,3,0,2,1,5,6,5,5,2,4,6,5,0,1,2,3,3,4,2,4,7,5)
cat(sprintf("   临床评分  mean: Ctrl=%.3f PRRSV=%.3f APRR=%.3f | SD: %.3f / %.3f / %.3f | n obs: %d/%d/%d\n",
            mean(ctrl_sc), mean(prr_sc), mean(apr_sc), sd(ctrl_sc), sd(prr_sc), sd(apr_sc),
            length(ctrl_sc), length(prr_sc), length(apr_sc)))
cat("   注意：以上 SD 为『每日观测值』的 SD，不是 6 头猪个体均值的 SD\n")
cat("   单猪层面（每猪 11 天评分）无法从扁平数组还原 -> 需原始台账\n\n")

cat("=== 4. 死亡 ===\n")
cat("   Control 0/6, PRRSV 2/6, APRR 1/6 (pzfx Table22 生存曲线)\n")
cat("   Fisher 精确检验 PRRSV vs APRR: ")
m <- matrix(c(2,4,1,5), nrow=2)
print(fisher.test(m)$p.value)
cat("\n=== 5. 体温 ===\n")
cat("   Control mean 39.467 (n=66 obs), PRRSV 40.672 (n=64), APRR 40.164 (n=66)\n")
cat("   -> 与临床评分同为重复测量\n")
