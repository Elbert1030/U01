#### ====================== 只展示，不输出 ====================== ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## 健壮性检查
seurat.vivo.cross<-readRDS("/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/PreProcessing/Harmony_integration_Cutoff_1K/cross.exp.combined.vivo.hto.celltag.annotation.rds")
md <- seurat.vivo.cross@meta.data
out_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1"
## 函数（保持不变）
get_clonal_matrix <- function(seurat, cluster_col, larry_col, sample_col,
                              samples_to_filt, min_cells, clusters_to_rm=NULL,
                              bcs_to_filt=NULL,
                              normalize_clusters=TRUE,
                              normalize_intraclone=TRUE,
                              normalize_pre_subset=FALSE) {
  larry_col_s   <- rlang::sym(larry_col)
  sample_col_s  <- rlang::sym(sample_col)
  cluster_col_s <- rlang::sym(cluster_col)
  
  clonal_mat <- FetchData(seurat, c(cluster_col, larry_col, sample_col)) %>% 
    dplyr::filter(!is.na(!!larry_col_s),
                  !!sample_col_s %in% samples_to_filt,
                  !(!!cluster_col_s) %in% clusters_to_rm) %>% 
    dplyr::select(-!!sample_col_s) %>% 
    dplyr::group_by(!!larry_col_s) %>%
    dplyr::filter(dplyr::n() >= min_cells) %>% 
    dplyr::ungroup() %>% 
    dplyr::count(!!cluster_col_s, !!larry_col_s) %>% 
    tidyr::pivot_wider(names_from = !!cluster_col_s, values_from = n, values_fill = 0) %>% 
    tibble::column_to_rownames(larry_col) %>% 
    as.matrix()
  
  # —— 安全子集：先子集后归一化（可选）——
  if (!is.null(bcs_to_filt) && normalize_pre_subset) {
    present <- bcs_to_filt[bcs_to_filt %in% rownames(clonal_mat)]
    clonal_mat <- clonal_mat[present, , drop = FALSE]
  }
  
  if (normalize_clusters) {
    clonal_mat <- apply(clonal_mat, 2, function(x) if (sum(x) == 0) x else x/sum(x))
  }
  
  if (normalize_intraclone) {
    clonal_mat <- t(apply(clonal_mat, 1, function(x) if (sum(x) == 0) x else x/sum(x)))
  }
  
  # —— 安全子集：先归一化后子集（默认）——
  if (!is.null(bcs_to_filt) && !normalize_pre_subset) {
    present <- bcs_to_filt[bcs_to_filt %in% rownames(clonal_mat)]
    clonal_mat <- clonal_mat[present, , drop = FALSE]
  }
  
  return(clonal_mat)
}




#### ===================================================== ####


## ========= 1) 找 OO / OY 的 clone 集合 & 共同克隆 =========
md <- seurat.vivo.cross@meta.data

oo_samples <- unique(md$orig.ident[grepl("OO", md$orig.ident)])
oy_samples <- unique(md$orig.ident[grepl("OY", md$orig.ident)])

# 各自样本中的唯一克隆
oo_clones_all <- unique(md$CloneID[
  !is.na(md$CloneID) & md$orig.ident %in% oo_samples
])
oy_clones_all <- unique(md$CloneID[
  !is.na(md$CloneID) & md$orig.ident %in% oy_samples
])

# 共同克隆
common_clones <- intersect(oo_clones_all, oy_clones_all)

# （默认）只保留 Old* 开头的共同克隆；如不需要此筛选，请注释下一行
common_clones <- common_clones[grepl("^Old", common_clones)]

length(common_clones); head(common_clones)

if (length(common_clones) == 0) stop("没有找到 OO 与 OY 的共同（Old*）克隆。")

## ========= 2) 构建两个矩阵：OO & OY（行=共同克隆；列=celltype） =========
min_cells <- 1  # 可调：提高到 2/3 能去噪

clonal_mat_OO <- get_clonal_matrix(
  seurat               = seurat.vivo.cross,
  cluster_col          = "celltype",
  larry_col            = "cloneTraceCombined",
  sample_col           = "orig.ident",
  samples_to_filt      = oo_samples,
  min_cells            = min_cells,
  bcs_to_filt          = common_clones,   # 只保留共同克隆
  normalize_clusters   = TRUE,            # 列归一化
  normalize_intraclone = TRUE,            # 行归一化
  normalize_pre_subset = TRUE             # 先子集（按共同克隆）再归一化（推荐）
)

clonal_mat_OY <- get_clonal_matrix(
  seurat               = seurat.vivo.cross,
  cluster_col          = "celltype",
  larry_col            = "cloneTraceCombined",
  sample_col           = "orig.ident",
  samples_to_filt      = oy_samples,
  min_cells            = min_cells,
  bcs_to_filt          = common_clones,
  normalize_clusters   = TRUE,
  normalize_intraclone = TRUE,
  normalize_pre_subset = TRUE
)

## ========= 3) 固定 x 轴（列）顺序；对齐两边的列和行 =========
desired_order <- c('LT-HSC', 'ST-HSC', 'MPP','MPP3', 'MPP4', 'GMP', 'MDP', 'MkP', 'MEP',
                   'EryP', 'CLP', 'DC', 'Mac', 'T cell', 'Granulocyte',
                   'UN1', 'UN2', 'UN3', 'UN4', 'UN5', 'UN6')

## ====== Keep only cell types that actually have data in OO or OY,
## ====== and show common OLD clone IDs as row labels on the heatmaps

## 1) Build union columns (desired_order first) and pad + reorder (you already did)
cols_union <- union(colnames(clonal_mat_OO), colnames(clonal_mat_OY))
cols_final <- c(intersect(desired_order, cols_union), setdiff(cols_union, desired_order))

pad_reorder <- function(M, cols) {
  miss <- setdiff(cols, colnames(M))
  if (length(miss)) {
    M <- cbind(M, matrix(0, nrow(M), length(miss),
                         dimnames = list(rownames(M), miss)))
  }
  M[, cols, drop = FALSE]
}
clonal_mat_OO <- pad_reorder(clonal_mat_OO, cols_final)
clonal_mat_OY <- pad_reorder(clonal_mat_OY, cols_final)

## 2) Row align = common clones (already “common Old* clones” from your pipeline)
rows_common <- intersect(rownames(clonal_mat_OO), rownames(clonal_mat_OY))
clonal_mat_OO <- clonal_mat_OO[rows_common, , drop = FALSE]
clonal_mat_OY <- clonal_mat_OY[rows_common, , drop = FALSE]

## 3) DROP columns (cell types) that are all-zero in BOTH OO and OY
keep_cols <- cols_final[colSums(clonal_mat_OO) > 0 | colSums(clonal_mat_OY) > 0]
clonal_mat_OO <- clonal_mat_OO[, keep_cols, drop = FALSE]
clonal_mat_OY <- clonal_mat_OY[, keep_cols, drop = FALSE]

## (optional) also drop clones (rows) that are all-zero in BOTH
keep_rows <- rowSums(clonal_mat_OO) > 0 | rowSums(clonal_mat_OY) > 0
clonal_mat_OO <- clonal_mat_OO[keep_rows, , drop = FALSE]
clonal_mat_OY <- clonal_mat_OY[keep_rows, , drop = FALSE]

## 4) Heatmaps: show clone IDs as row labels; fix column_order to keep_cols
vmax <- as.numeric(quantile(c(clonal_mat_OO, clonal_mat_OY), 0.98, na.rm = TRUE))
col_fun <- circlize::colorRamp2(c(0, vmax), hcl_palette = "Rocket", reverse = TRUE)

ht_OO <- Heatmap(
  clonal_mat_OO,
  name = "fraction",
  col  = col_fun,
  show_row_names   = FALSE,                     # ← show clone IDs
  row_names_gp     = gpar(fontsize = 10),
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = TRUE,                     # OO determines row order
  cluster_columns  = FALSE,
  column_order     = keep_cols,                # ← only non-empty cell types
  border           = TRUE,
  width  = unit(4.5, "inch"),
  height = unit(4.2, "inch"),
  column_title = "OO samples (common Old* clones)",
  row_title     = paste0("Common Old* clones (n=", nrow(clonal_mat_OO), ")"),
  column_names_gp = gpar(fontsize = 11)
)

ht_OY <- Heatmap(
  clonal_mat_OY,
  name = "fraction",
  col  = col_fun,
  show_row_names   = FALSE,                     # ← show clone IDs
  row_names_gp     = gpar(fontsize = 10),
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,                    # follow OO’s row order
  cluster_columns  = FALSE,
  column_order     = keep_cols,                # same kept columns/order
  border           = TRUE,
  width  = unit(4.5, "inch"),
  height = unit(4.2, "inch"),
  column_title = "OY samples (common Old* clones)",
  column_names_gp = gpar(fontsize = 11)
)

draw(ht_OO + ht_OY, ht_gap = unit(6, "mm"), main_heatmap = 1)

## === OO: CloneID × (Rep1 / Rep2 / Rep3 / UnMapped) ===
## 结果同时给“计数矩阵”和“行内比例矩阵”，行顺序与 clonal_mat_OO 一致

# 1) 取要展示的克隆（与热图一致）
clones_oo <- rownames(clonal_mat_OO)

# 2) 从 meta 中抽取 OO 的细胞并规整 RepExtended
md_oo <- md[md$orig.ident %in% oo_samples &
              !is.na(md$CloneID) &
              md$CloneID %in% clones_oo,
            c("CloneID","Rep")]
md_oo$Rep <- trimws(as.character(md_oo$Rep))

# 3) 归类为四列：Rep1/Rep2/Rep3/UnMapped（只要不是 OO_Rep1/2/3 就归到 UnMapped）
md_oo$rep_bin <- dplyr::case_when(
  md_oo$Rep == "Rep1" ~ "Rep1",
  md_oo$Rep == "Rep2" ~ "Rep2",
  md_oo$Rep == "Rep3" ~ "Rep3",
  TRUE                            ~ "UnMapped"
)

# 4) 计数矩阵（缺失列补 0），并按热图的克隆顺序排列
tbl_counts_OO <- md_oo |>
  dplyr::count(CloneID, rep_bin, name = "cells") |>
  tidyr::pivot_wider(names_from = rep_bin, values_from = cells, values_fill = 0) |>
  tibble::column_to_rownames("CloneID")

need_cols <- c("Rep1","Rep2","Rep3","UnMapped")
miss <- setdiff(need_cols, colnames(tbl_counts_OO))
if (length(miss)) {
  tbl_counts_OO[, miss] <- 0
}
tbl_counts_OO <- tbl_counts_OO[clones_oo, need_cols, drop = FALSE]

# 5) 行内比例矩阵（每个克隆在四列上的分布）
rs <- rowSums(tbl_counts_OO)
tbl_frac_OO <- sweep(tbl_counts_OO, 1, ifelse(rs==0, 1, rs), "/")

# 6) 打印（需要可保存为 CSV）
cat("\n== OO: Clone × Rep counts ==\n")
print(tibble::as_tibble(cbind(clone_id = rownames(tbl_counts_OO),
                              as.data.frame(tbl_counts_OO))), n = Inf)

cat("\n== OO: Clone × Rep fractions (row-normalized) ==\n")
print(tibble::as_tibble(cbind(clone_id = rownames(tbl_frac_OO),
                              as.data.frame(round(tbl_frac_OO, 4)))), n = Inf)

# 可选：保存
 write.csv(cbind(clone_id = rownames(tbl_counts_OO), as.data.frame(tbl_counts_OO)),
           file = file.path(out_dir,"OO_clone_by_rep_counts.csv"), row.names = FALSE)
 write.csv(cbind(clone_id = rownames(tbl_frac_OO), as.data.frame(tbl_frac_OO)),
           file = file.path(out_dir,"OO_clone_by_rep_fractions.csv"), row.names = FALSE)
## === OY: CloneID × (Rep1 / Rep2 / Rep3 / UnMapped) ===
## 与热图一致的行顺序；输出计数矩阵与行内比例矩阵

# 1) 克隆顺序（与热图一致）
clones_oy <- rownames(clonal_mat_OY)

# 2) 取 OY 细胞并规整 RepExtended
md_oy <- md[md$orig.ident %in% oy_samples &
              !is.na(md$CloneID) &
              md$CloneID %in% clones_oy,
            c("CloneID","Rep")]
md_oy$Rep <- trimws(as.character(md_oy$Rep))

# 3) 映射到四列：Rep1/Rep2/Rep3/UnMapped
md_oy$rep_bin <- dplyr::case_when(
  md_oy$Rep == "Rep1" ~ "Rep1",
  md_oy$Rep == "Rep2" ~ "Rep2",
  md_oy$Rep == "Rep3" ~ "Rep3",
  TRUE                            ~ "UnMapped"
)

# 4) 计数矩阵（缺失列补 0），并按克隆顺序排列
tbl_counts_OY <- md_oy |>
  dplyr::count(CloneID, rep_bin, name = "cells") |>
  tidyr::pivot_wider(names_from = rep_bin, values_from = cells, values_fill = 0) |>
  tibble::column_to_rownames("CloneID")

need_cols <- c("Rep1","Rep2","Rep3","UnMapped")
miss <- setdiff(need_cols, colnames(tbl_counts_OY))
if (length(miss)) tbl_counts_OY[, miss] <- 0
tbl_counts_OY <- tbl_counts_OY[clones_oy, need_cols, drop = FALSE]

# 5) 行内比例矩阵
rs <- rowSums(tbl_counts_OY)
tbl_frac_OY <- sweep(tbl_counts_OY, 1, ifelse(rs==0, 1, rs), "/")

# 6) 打印（需要可保存为 CSV）
cat("\n== OY: Clone × Rep counts ==\n")
print(tibble::as_tibble(cbind(clone_id = rownames(tbl_counts_OY),
                              as.data.frame(tbl_counts_OY))), n = Inf)

cat("\n== OY: Clone × Rep fractions (row-normalized) ==\n")
print(tibble::as_tibble(cbind(clone_id = rownames(tbl_frac_OY),
                              as.data.frame(round(tbl_frac_OY, 4)))), n = Inf)

# 可选保存：
 write.csv(cbind(clone_id = rownames(tbl_counts_OY), as.data.frame(tbl_counts_OY)),
           file = file.path(out_dir,"OY_clone_by_rep_counts.csv"), row.names = FALSE)
 write.csv(cbind(clone_id = rownames(tbl_frac_OY), as.data.frame(tbl_frac_OY)),
           file = file.path(out_dir,"OY_clone_by_rep_fractions.csv"), row.names = FALSE)

library(ComplexHeatmap)

## === Step 1: 主要 replicate (每行最大值) ===
rep_order <- c("Rep1","Rep2","Rep3","UnMapped")
row_rep_OO <- apply(tbl_frac_OO, 1, function(x) rep_order[which.max(x)])
row_rep_OO <- factor(row_rep_OO, levels = rep_order)

## === Step 2: 定义颜色 ===
rep_colors <- c(
  "Rep1"     = "#1f77b4",   # 蓝
  "Rep2"     = "#2ca02c",   # 绿
  "Rep3"     = "#ff7f0e",   # 橙
  "UnMapped" = "#7f7f7f"    # 灰
)

## === Step 3: 按照 replicate 顺序排序行 ===
row_order <- order(row_rep_OO)   # factor 保证 Rep1→Rep2→Rep3→UnMapped
clonal_mat_OO <- clonal_mat_OO[row_order, , drop = FALSE]
row_rep_OO    <- row_rep_OO[row_order]

## === Step 4: 构建行注释 ===
row_ha_OO <- rowAnnotation(
  Replicate = row_rep_OO,
  col = list(Replicate = rep_colors),
  width = unit(0.6, "cm")
)

## === Step 5: Heatmap 本体 ===
ht_OO <- Heatmap(
  clonal_mat_OO,
  name = "fraction",
  col  = col_fun,
  show_row_names   = FALSE,
  row_names_gp     = gpar(fontsize = 6),
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,                  # ⚠️ 这里要 FALSE，否则聚类会打乱顺序
  cluster_columns  = FALSE,
  column_order     = keep_cols,
  border           = TRUE,
  width  = unit(4.5, "inch"),
  height = unit(4.2, "inch"),
  column_title = "OO samples (common Old* clones)",
  row_title     = paste0("Common Old* clones (n=", nrow(clonal_mat_OO), ")"),
  column_names_gp = gpar(fontsize = 11)
)

## === Step 6: 绘图 ===
draw(row_ha_OO + ht_OO,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")


library(ComplexHeatmap)

## === Step 1: 主要 replicate (每行最大值) ===
rep_order <- c("Rep1","Rep2","Rep3","UnMapped")
row_rep_OY <- apply(tbl_frac_OY, 1, function(x) rep_order[which.max(x)])
row_rep_OY <- factor(row_rep_OY, levels = rep_order)

## === Step 2: 定义颜色 ===
rep_colors <- c(
  "Rep1"     = "#1f77b4",   # 蓝
  "Rep2"     = "#2ca02c",   # 绿
  "Rep3"     = "#ff7f0e",   # 橙
  "UnMapped" = "#7f7f7f"    # 灰
)

## === Step 3: 按照 replicate 顺序排序行 ===
row_order_OY <- order(row_rep_OY)   # factor 确保 Rep1→Rep2→Rep3→UnMapped
clonal_mat_OY <- clonal_mat_OY[row_order_OY, , drop = FALSE]
row_rep_OY    <- row_rep_OY[row_order_OY]

## === Step 4: 构建行注释 ===
row_ha_OY <- rowAnnotation(
  Replicate = row_rep_OY,
  col = list(Replicate = rep_colors),
  width = unit(0.6, "cm")
)

## === Step 5: Heatmap 本体 ===
ht_OY <- Heatmap(
  clonal_mat_OY,
  name = "fraction",
  col  = col_fun,
  show_row_names   = FALSE,
  row_names_gp     = gpar(fontsize = 8),
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,    # ⚠️ 防止聚类打乱 Rep 顺序
  cluster_columns  = FALSE,
  column_order     = keep_cols,
  border           = TRUE,
  column_title     = "OY samples (common Old* clones)"
)

## === Step 6: 绘图 ===
draw(row_ha_OY + ht_OY,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

## === Step 1: 定义顺序（基于 YY 的行次序） ===
# row_rep_YY 是之前算出来的 replicate 因子
# 用它的排序作为 "模板"
row_order_template <- order(row_rep_OO)

# 按这个顺序重排 YY
clonal_mat_OO <- clonal_mat_OO[row_order_template, , drop = FALSE]
row_rep_OO    <- row_rep_OO[row_order_template]

# 按相同的行名顺序重排 YO（确保 cloneID 对齐）
clonal_mat_OY <- clonal_mat_OY[rownames(clonal_mat_OO), , drop = FALSE]

## === Step 2: YO replicate 信息（用 YO 的表重新计算） ===
rep_order <- c("Rep1","Rep2","Rep3","UnMapped")
row_rep_OY <- apply(tbl_frac_OY, 1, function(x) rep_order[which.max(x)])
row_rep_OY <- factor(row_rep_OY, levels = rep_order)
row_rep_OY <- row_rep_OY[rownames(clonal_mat_OY)]   # 保持和行顺序一致

## === Step 3: 行注释 ===
rep_colors <- c("Rep1"="#1f77b4","Rep2"="#2ca02c","Rep3"="#ff7f0e","UnMapped"="#7f7f7f")

row_ha_OO <- rowAnnotation(
  Replicate = row_rep_OO,
  col = list(Replicate = rep_colors),
  width = unit(0.6, "cm")
)

row_ha_OY <- rowAnnotation(
  Replicate = row_rep_OY,
  col = list(Replicate = rep_colors),
  width = unit(0.6, "cm")
)

## === Step 4: Heatmaps ===
ht_OO <- Heatmap(
  clonal_mat_OO, name="fraction", col=col_fun,
  show_row_names=FALSE, row_names_gp=gpar(fontsize=8),
  cluster_rows=FALSE, cluster_columns=FALSE,
  column_order=keep_cols, border=TRUE,
  column_title="OO samples (common Old* clones)",
  row_title=paste0("Common Old* clones (n=", nrow(clonal_mat_OO), ")")
)

ht_OY <- Heatmap(
  clonal_mat_OY, name="fraction", col=col_fun,
  show_row_names=FALSE, row_names_gp=gpar(fontsize=8),
  cluster_rows=FALSE, cluster_columns=FALSE,
  column_order=keep_cols, border=TRUE,
  column_title="OY samples (common Old* clones)"
)

## === Step 5: 并排画图 ===
ht_all <- (row_ha_OO + ht_OO) + (row_ha_OY + ht_OY)

draw(ht_all,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

# 定义输出目录
#out_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 输出 PNG
#png(file.path(out_dir, "OO_OY_heatmap.png"), width = 12, height = 6, units = "in", res = 300)
#draw((row_ha_OO + ht_OO) + (row_ha_OY + ht_OY),
#     heatmap_legend_side = "right",
#     annotation_legend_side = "right")
#dev.off()

# 输出 PDF
#pdf(file.path(out_dir, "OO_OY_heatmap.pdf"), width = 12, height = 6)
#draw((row_ha_OO + ht_OO) + (row_ha_OY + ht_OY),
#     heatmap_legend_side = "right",
#     annotation_legend_side = "right")
#dev.off()
library(dplyr)
library(ggplot2)

## 1) 用 clone × replicate 的 count 矩阵求总 clone size
clone_size_OO <- rowSums(tbl_counts_OO)
clone_size_OY <- rowSums(tbl_counts_OY)

## 2) 合并成 dataframe
df <- data.frame(
  clone_id = intersect(names(clone_size_OO), names(clone_size_OY)),
  size_OO  = clone_size_OO[intersect(names(clone_size_OO), names(clone_size_OY))],
  size_OY  = clone_size_OY[intersect(names(clone_size_OO), names(clone_size_OY))]
)

## 3) 计算 log2FC (OO vs OY)，加 pseudocount=1 防止 log2(0)
df <- df %>%
  mutate(
    log2FC    = log2((size_OO + 1) / (size_OY + 1)),
    cloneSize = (size_OO + size_OY) / sum(size_OO + size_OY)  # 占比用于点大小
  )

## 4) 按 log2FC 排序并加 rank
df <- df %>%
  arrange(log2FC) %>%
  mutate(rank = row_number())

## 5) 绘图
p <- ggplot(df, aes(x = rank, y = log2FC)) +
  geom_point(aes(size = cloneSize, color = log2FC)) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0,
                        limits=c(-6,6)) +
  scale_size(range = c(1, 8)) +
  theme_classic(base_size = 16) +
  labs(
    title = "OO vs OY Clone Size Comparison",
    x = "Rank",
    y = "log2FC (OO vs OY)",
    color = "log2FC",
    size = "Clone size"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-2.5, 2.5), linetype = "dashed") +
  annotate("text", x = max(df$rank)*0.1, y = 3, label = "Clone size increased", color = "red") +
  annotate("text", x = max(df$rank)*0.1, y = -3, label = "Clone size decreased", color = "blue")

print(p)

# 保存结果
#out_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
#ggsave(file.path(out_dir, "OO_vs_OY_cloneSize_rankedLog2FC.png"), p, width=7, height=6, dpi=300)

