# =============================================================================
# CellChat LR network analysis — reproduce old Fig6 (APRR anti-PRRSV, 3-group)
#   Control / PRRSV / APRR4
# Workflow (CellChat 1.5.0):
#   createCellChat -> subsetData -> identifyOverExpressedGenes
#     -> identifyOverExpressedInteractions  (BUILD @LR$LRsig — REQUIRED)
#     -> computeCommunProb -> filterCommunication
#     -> computeCommunProbPathway -> aggregateNet -> viz
# =============================================================================
suppressMessages(library(Seurat))
suppressMessages(library(CellChat))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
# 顺序执行：computeCommunProb 的 bootstrap 在主进程运行，确保我们的 namespace 覆盖生效
plan(sequential)
set.seed(42)   # reproducibility for computeCommunProb bootstrap
cat(">>> future plan:", class(future::plan()), "\n")

# ---------------------------------------------------------------------------
# 不再覆盖 CellChat 内部 computeExpr_* 函数。
# 原因：之前用 dplyr 的兜底 guard 在热循环里反复调用 dplyr，已知会触发
# Rscript+dplyr 的 segfault。现在 DB 已清洗(所有复合物亚基都在 data 里)、
# @LR$LRsig 已正确构建(775 对、过滤丢 0 条、配/受体全在 data)，原始函数可直接跑。
# 仅保留一个 base-R 的轻微兜底，避免任何 NA 亚基导致的 geometricMean 崩溃。
# ---------------------------------------------------------------------------
.computeExpr_complex_safe <- function(complex_input, data.use, complex) {
  sub_cols <- grep("^subunit", colnames(complex_input), value = TRUE)
  RV <- unlist(as.list(complex_input[complex, sub_cols, drop = FALSE]), use.names = FALSE)
  RV <- RV[RV != "" & !is.na(RV)]
  RV <- RV[RV %in% rownames(data.use)]
  if (length(RV) == 0) return(matrix(NA_real_, nrow = length(complex), ncol = ncol(data.use)))
  CellChat:::geometricMean(data.use[RV, , drop = FALSE])
}
assignInNamespace("computeExpr_complex", .computeExpr_complex_safe, ns = "CellChat")
# 不覆盖 computeExpr_LR（保留 CellChat 原版）；仅 computeExpr_complex 用 base-R 兜底以防 NA 亚基。
cat(">>> computeExpr_complex base-R guard active (no dplyr)\n")

# flush helper for visible progress
.flog <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

# ---------------------------------------------------------------------------
# igraph 2.x patch: i.parse.plot.params errors with "replacement has length zero"
# when a vertex/edge attribute has NA values but igraph 2.x has NO registered default
# (i.default.values[[type]][[name]] is NULL). CellChat 1.5.0 triggers this in
# netVisual_circle (and any igraph-based CellChat plot) under igraph >= 2.0.
# We patch i.parse.plot.params so a NULL default is replaced with NA (a legal no-op)
# instead of crashing. This is a runtime-only patch; the shared R library is untouched.
# ---------------------------------------------------------------------------
# Safe fallback when igraph 2.x has no registered default for an attribute:
# numeric-ish attributes -> 1, everything else -> "" (avoids NA in logical tests
# inside plot.igraph which would throw "missing value where TRUE/FALSE needed").
.igraph_dv_fallback <- function(name) {
  # Any *label* attribute must default to NULL: plot.igraph guards these with
  # `!is.null(...)`, so NULL is skipped safely. Returning a string ("black") would
  # be fed into `if(attr)` and throw "argument is not interpretable as logical".
  if (grepl("label", name, ignore.case = TRUE)) return(NULL)         # label* -> NULL
  if (name == "frame") return(FALSE)                                  # plot.frame is logical
  if (grepl("(color|colour|fill|lty|lwd|arrow)", name, ignore.case = TRUE)) return("black")  # color attrs
  if (grepl("(size|width|dist|angle|adj|font|cex|curv|loop|scaling|shape|margin|layout|rescale)", name, ignore.case = TRUE)) return(1)  # numeric attrs
  return("black")                                                     # safe string fallback
}
.i.parse.plot.params_patched <- function(graph, params) {
  p <- list(vertex = list(), edge = list(), plot = list())
  for (n in names(params)) {
    if (substr(n, 1, 7) == "vertex.") {
      nn <- substring(n, 8); p[["vertex"]][[nn]] <- params[[n]]
    } else if (substr(n, 1, 5) == "edge.") {
      nn <- substring(n, 6); p[["edge"]][[nn]] <- params[[n]]
    } else {
      p[["plot"]][[n]] <- params[[n]]
    }
  }
  # igraph 2.x regression: plotting flags that reach plot.default as logical
  # (frame/frame.plot/axes/rescale) can arrive as non-logical (string/NA) and
  # crash `if (frame.plot)`. Force them logical up front.
  for (ln in c("frame", "frame.plot", "axes", "rescale")) {
    if (!is.null(p[["plot"]][[ln]])) {
      if (!is.logical(p[["plot"]][[ln]]) || any(is.na(p[["plot"]][[ln]]))) {
        p[["plot"]][[ln]] <- FALSE
      }
    }
  }
  func <- function(type = c("vertex", "edge", "plot"), name, range = NULL, dontcall = FALSE) {
    type <- igraph:::igraph_match_arg(type)
    # Intercept logical plot flags BEFORE the direct-return bypass below.
    if (type == "plot" && (name %in% c("frame", "frame.plot", "axes", "rescale"))) {
      v <- p[["plot"]][[name]]
      if (is.logical(v) && !any(is.na(v))) return(v)
      return(FALSE)
    }
    ret <- function() {
      v <- p[[type]][[name]]
      if (is.function(v) && !dontcall) v <- v(graph)
      if (is.null(range)) return(v)
      else {
        if (length(v) == 1) return(rep(v, length(range)))
        else return(rep(v, length.out = max(range) + 1)[[range + 1]])
      }
    }
    if (name %in% names(p[[type]])) return(ret())
    else {
      if (type == "vertex" && name %in% igraph:::vertex_attr_names(graph)) {
        p[[type]][[name]] <- igraph:::vertex_attr(graph, name)
      } else if (type == "edge" && name %in% igraph:::edge_attr_names(graph)) {
        p[[type]][[name]] <- igraph:::edge_attr(graph, name)
      } else if (type == "plot" && name %in% igraph:::graph_attr_names(graph)) {
        p[[type]][[name]] <- igraph:::graph_attr(graph, name)
      } else {
        n <- paste(sep = "", type, ".", name)
        v <- igraph:::igraph_opt(n)
        if (!is.null(v)) { p[[type]][[name]] <- v; return(ret()) }
        dv <- igraph:::i.default.values[[type]][[name]]
        if (is.null(dv)) dv <- .igraph_dv_fallback(name)
        p[[type]][[name]] <- dv
        return(ret())
      }
    }
    if (!is.function(p[[type]][[name]])) {
      if (any(is.na(p[[type]][[name]]))) {
        if (name != "label") {
          dv <- igraph:::i.default.values[[type]][[name]]
          if (is.null(dv)) dv <- .igraph_dv_fallback(name)
          p[[type]][[name]][is.na(p[[type]][[name]])] <- dv
        } else {
          p[[type]][[name]][is.na(p[[type]][[name]])] <- ""
        }
      }
    }
    return(ret())
  }
  return(func)
}
assignInNamespace("i.parse.plot.params", .i.parse.plot.params_patched, ns = "igraph")
cat(">>> igraph i.parse.plot.params patched (igraph 2.x compat)\n")

setwd("data/scRNA-seq/reanalysis_v2")
cat(">>> loading 18ct rds ...\n")
obj <- readRDS("data/derived/seurat_reanalysis_v2_18ct.rds")
obj <- JoinLayers(obj, assay = "RNA")
# NOTE: this Seurat object stores ONLY raw counts (the 'data' layer is empty after
# JoinLayers). The original Fig6 was therefore built from raw counts. CellChat accepts
# the raw count matrix as @data (it scales internally via data/max(data)), so we pass
# counts directly to reproduce the old figures exactly.
counts <- LayerData(obj, assay = "RNA", layer = "counts")
data.input <- as(counts, "dgCMatrix")
cat("data.input dims:", nrow(data.input), "x", ncol(data.input), "\n")
gn <- rownames(data.input)

# ---- Clean the GLOBAL CellChatDB.human (subsetData rebuilds @DB from it) ----
cpx <- CellChatDB.human$complex
cpx_cols <- paste0("subunit_", 1:4)
cpx_names <- rownames(cpx)
cpx_keep <- apply(as.matrix(cpx[, cpx_cols, drop = FALSE]), 1,
                  function(r) all(r[r != "" & !is.na(r)] %in% gn))
CellChatDB.human$complex <- cpx[cpx_keep, , drop = FALSE]
cat("complex kept:", nrow(CellChatDB.human$complex), "of", nrow(cpx), "\n")

keep_int <- !is.na(CellChatDB.human$interaction$ligand) & !is.na(CellChatDB.human$interaction$receptor) &
  (CellChatDB.human$interaction$ligand %in% gn | CellChatDB.human$interaction$ligand %in% cpx_names) &
  (CellChatDB.human$interaction$receptor %in% gn | CellChatDB.human$interaction$receptor %in% cpx_names)
for (gc in c("agonist", "antagonist")) {
  if (gc %in% colnames(CellChatDB.human$interaction)) {
    v <- CellChatDB.human$interaction[[gc]]
    keep_int <- keep_int & (is.na(v) | v == "" | v %in% gn)
  }
}
cat("interactions kept:", sum(keep_int), "of", nrow(CellChatDB.human$interaction), "\n")
CellChatDB.human$interaction <- CellChatDB.human$interaction[keep_int, , drop = FALSE]

cf <- CellChatDB.human$cofactor
cf_cols <- paste0("cofactor", 1:16)
cf_ok <- apply(as.matrix(cf[, cf_cols, drop = FALSE]), 1, function(r) all(r[r != ""] %in% gn))
CellChatDB.human$cofactor <- cf[cf_ok, , drop = FALSE]
cat("cofactor kept:", nrow(CellChatDB.human$cofactor), "\n")

make_db <- function() {
  d <- CellChatDB.human
  d$interaction <- data.frame(d$interaction)
  d$complex <- data.frame(d$complex)
  d$cofactor <- data.frame(d$cofactor)
  d
}
db <- make_db()

# --- SPEED FIX: subset input to ONLY genes referenced by the (cleaned) DB ---
# computeCommunProb internally does `as.matrix(object@data.signaling)` -> a DENSE matrix,
# then `t()` + 100 bootstrap `aggregate()` calls over ALL genes x cells (~3.5 GB, extremely
# slow / appears hung). CellChat only EVER indexes DB genes for LR computation, so non-DB
# genes contribute NOTHING to the result. Subsetting to DB genes shrinks the dense matrix
# ~9x with IDENTICAL LR probabilities.
sub_cols <- grep("^subunit", colnames(db$complex), value = TRUE)
dbg <- c(as.character(db$interaction$ligand), as.character(db$interaction$receptor),
         as.character(db$interaction$agonist), as.character(db$interaction$antagonist),
         unlist(as.list(db$complex[, sub_cols, drop = FALSE]), use.names = FALSE))
dbg <- dbg[!is.na(dbg) & dbg != "" & dbg %in% gn]
data.input <- data.input[rownames(data.input) %in% dbg, , drop = FALSE]
gn <- rownames(data.input)
cat(">>> data.input subset to DB genes:", nrow(data.input), "x", ncol(data.input), "\n")

# fix_signaling: ensure @data.signaling / @data.project hold ALL (subset) genes so any
# DB-referenced gene (incl. complex subunits / cofactors) is present for computeCommunProb.
fix_signaling <- function(cc) {
  cc@data.signaling <- cc@data
  cc@data.project  <- cc@data
  cc
}

meta <- data.frame(
  cell_type = obj$cell_type,
  group     = obj$group,
  row.names = colnames(obj)
)

cat("=== building CellChat object (group) ===\n")
cellchat <- createCellChat(data.input, meta = meta, group.by = "group")
cellchat@DB <- db
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)   # <-- BUILDS @LR$LRsig
cellchat <- fix_signaling(cellchat)
cellchat@DB <- make_db()
cat(">>> @LR$LRsig built: nrow =", nrow(cellchat@LR$LRsig), "\n")

# safety filter: drop NA or gene/complex-name missing from data
lrk <- !is.na(cellchat@LR$LRsig$ligand) & !is.na(cellchat@LR$LRsig$receptor) &
  (cellchat@LR$LRsig$ligand %in% gn | cellchat@LR$LRsig$ligand %in% cpx_names) &
  (cellchat@LR$LRsig$receptor %in% gn | cellchat@LR$LRsig$receptor %in% cpx_names)
cat(">>> @LR$LRsig filter: dropped", sum(!lrk), "of", nrow(cellchat@LR$LRsig), "\n")
cellchat@LR$LRsig <- cellchat@LR$LRsig[lrk, ]

.flog(">>> ENTER computeCommunProb (group) ...\n")
cellchat <- tryCatch(
  computeCommunProb(cellchat, type = "triMean", trim = 0.1, do.fast = TRUE),
  error = function(e) { .flog("!!! computeCommunProb(group) ERROR: %s\n", conditionMessage(e)); stop(e) })
.flog(">>> computeCommunProb (group) OK\n")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cat(">>> group CellChat net built; pathways =", length(cellchat@netP$pathways), "\n")
saveRDS(cellchat, "data/derived/cellchat_group.rds")

W1 <- 3.5
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(4, 4, 4, 4, "mm"))
}

# Fig6A: global circle plot (3-group)
tryCatch({
  pdf("_archive_legacy/figures_paper_3grp_panels/Fig6A_global_comm_circle.pdf", width = W1 * 1.6, height = W1 * 1.6)
  netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                   weight.scale = TRUE, label.edge = FALSE, title.name = "CellChat global — group")
  dev.off()
  png("_archive_legacy/figures_paper_3grp_panels/Fig6A_global_comm_circle.png", width = W1 * 1.6, height = W1 * 1.6, units = "in", res = 300)
  netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                   weight.scale = TRUE, label.edge = FALSE, title.name = "CellChat global — group")
  dev.off()
  cat(">>> Fig6A done\n")
}, error = function(e) cat("!!! Fig6A FAILED:", conditionMessage(e), "\n"))

# Fig6B: chord plot (global, 3-group). aggregateNet collapses LR pairs into a 3x3
# count/weight matrix stored in @net$count; pass it as `net` to the chord visualizer.
tryCatch({
  pdf("_archive_legacy/figures_paper_3grp_panels/Fig6B_global_chord.pdf", width = W1 * 1.6, height = W1 * 1.6)
  netVisual_chord_cell(cellchat, net = cellchat@net$count, title.name = "CellChat chord — group")
  dev.off()
  cat(">>> Fig6B done\n")
}, error = function(e) cat("!!! Fig6B FAILED:", conditionMessage(e), "\n"))

# Fig6C: heatmap of signaling pathways
tryCatch({
  pdf("_archive_legacy/figures_paper_3grp_panels/Fig6C_pathway_heatmap.pdf", width = W1 * 1.6, height = W1 * 1.8)
  netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "YlOrRd")
  dev.off()
  png("_archive_legacy/figures_paper_3grp_panels/Fig6C_pathway_heatmap.png", width = W1 * 1.6, height = W1 * 1.8, units = "in", res = 300)
  netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "YlOrRd")
  dev.off()
  cat(">>> Fig6C done\n")
}, error = function(e) cat("!!! Fig6C FAILED:", conditionMessage(e), "\n"))

# Focus key LR pathways: PD-L1/PD-1, SPP1, MIF
focus_lr <- c("PD-L1/PD-1", "SPP1", "MIF")
avail_path <- unique(cellchat@LR$LRsig$pathway_name)
for (lr in focus_lr) {
  if (lr %in% avail_path) {
    base <- paste0("_archive_legacy/figures_paper_3grp_panels/Fig6_", gsub("/", "_", lr), "_lr")
    tryCatch({
      pdf(paste0(base, ".pdf"), width = W1 * 1.4, height = W1 * 1.4)
      netVisual_individual(cellchat, signaling = lr, layout = "circle")
      dev.off()
      png(paste0(base, ".png"), width = W1 * 1.4, height = W1 * 1.4, units = "in", res = 300)
      netVisual_individual(cellchat, signaling = lr, layout = "circle")
      dev.off()
      cat("focus LR plotted:", lr, "\n")
    }, error = function(e) cat("focus LR FAILED:", lr, "->", conditionMessage(e), "\n"))
  } else {
    cat("LR pathway not in network:", lr, "\n")
  }
}

# ---- Cell-type resolved LR network (neutrophil vs T-cell PD-L1/PD-1 etc.) ----
cat("=== building CellChat object (cell_type) ===\n")
cellchat_ct <- createCellChat(data.input, meta = meta, group.by = "cell_type")
cellchat_ct@DB <- db
cellchat_ct <- subsetData(cellchat_ct)
cellchat_ct <- identifyOverExpressedGenes(cellchat_ct)
cellchat_ct <- identifyOverExpressedInteractions(cellchat_ct)   # <-- BUILDS @LR$LRsig
cellchat_ct <- fix_signaling(cellchat_ct)
cellchat_ct@DB <- make_db()
lrk_ct <- !is.na(cellchat_ct@LR$LRsig$ligand) & !is.na(cellchat_ct@LR$LRsig$receptor) &
  (cellchat_ct@LR$LRsig$ligand %in% gn | cellchat_ct@LR$LRsig$ligand %in% cpx_names) &
  (cellchat_ct@LR$LRsig$receptor %in% gn | cellchat_ct@LR$LRsig$receptor %in% cpx_names)
cellchat_ct@LR$LRsig <- cellchat_ct@LR$LRsig[lrk_ct, ]
cat(">>> celltype @LR$LRsig built: nrow =", nrow(cellchat_ct@LR$LRsig), "\n")
.flog(">>> ENTER computeCommunProb (cell_type) ...\n")
cellchat_ct <- tryCatch(
  computeCommunProb(cellchat_ct, type = "triMean", trim = 0.1, do.fast = TRUE),
  error = function(e) { .flog("!!! computeCommunProb(cell_type) ERROR: %s\n", conditionMessage(e)); stop(e) })
.flog(">>> computeCommunProb (cell_type) OK\n")
cellchat_ct <- filterCommunication(cellchat_ct, min.cells = 10)
cellchat_ct <- aggregateNet(cellchat_ct)
saveRDS(cellchat_ct, "data/derived/cellchat_ct.rds")
tryCatch({
  pdf("_archive_legacy/figures_paper_3grp_panels/Fig6D_celltype_interaction.pdf", width = W1 * 1.8, height = W1 * 1.8)
  netVisual_heatmap(cellchat_ct, measure = "weight", color.heatmap = "YlOrRd")
  dev.off()
  png("_archive_legacy/figures_paper_3grp_panels/Fig6D_celltype_interaction.png", width = W1 * 1.8, height = W1 * 1.8, units = "in", res = 300)
  netVisual_heatmap(cellchat_ct, measure = "weight", color.heatmap = "YlOrRd")
  dev.off()
  cat(">>> Fig6D done\n")
}, error = function(e) cat("!!! Fig6D FAILED:", conditionMessage(e), "\n"))

cat(">>> Fig6 CellChat DONE\n")
