#A Structure-Informed Atlas of Venom-Derived Peptides Reveals the Organization of Chemical Space
#Thaís Caroline Gonçalves, Eduardo Henrique Toral Cortez, and Danilo Trabuco Amaral

suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(ggrepel); library(ggpubr)
  library(ggcorrplot); library(uwot); library(dbscan); library(vegan)
  library(rstatix); library(ComplexHeatmap); library(circlize); library(broom)
  library(stringr); library(readr); library(MASS); library(randomForest)
})
set.seed(123)

# Optional: set a working directory if needed
setwd("C:/Users/Pichau/Desktop/Thais_C/")

in_tsv  <- "attrubutes_prot.txt"
outroot <- "atlas_out"
dir.create(outroot, showWarnings = FALSE)
dir.create(file.path(outroot, "figs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outroot, "tables"), showWarnings = FALSE)
dir.create(file.path(outroot, "supplementary"), showWarnings = FALSE)

# ==== 1) Input and cleaning ====
df0 <- readr::read_tsv(in_tsv, show_col_types = FALSE) |> janitor::clean_names()

df <- df0 |>
  mutate(
    comprimento     = coalesce(length, n_residues),
    # 'charge_pH7' becomes 'charge_p_h7' after clean_names()
    carga           = coalesce(!!rlang::sym("charge_p_h7"), charge_p_h7),
    frac_ss_helix   = coalesce(sec_helix, ss_pred_helix),
    frac_ss_sheet   = coalesce(sec_sheet, ss_pred_sheet),
    frac_ss_turn    = ss_pred_turn,
    classe = case_when(
      str_detect(input_file, regex("conotoxin|alpha_conotoxin|conorfamide", TRUE)) ~ "Conotoxin-like",
      str_detect(input_file, regex("phospholipase", TRUE))                          ~ "Phospholipase A2-like",
      str_detect(input_file, regex("bradykinin|kinin", TRUE))                       ~ "Bradykinin-related",
      str_detect(input_file, regex("neurotoxin|toxin", TRUE))                       ~ "Neurotoxin/Tx",
      str_detect(input_file, regex("antimicrobial", TRUE))                          ~ "AMP",
      TRUE ~ "Others"
    ),
    taxon_guess = str_extract(input_file, "(Conus|Conasprella|Bothrops|Daboia|Micrurus|Phoneutria|Physalaemus|Heterometrus|Lachesis|Tityus|Deinagkistrodon)"),
    n_cys  = cys_count,
    n_ss   = disulfide_pairs,
    has_signal = if_else(str_to_upper(signalp) == "YES", TRUE, FALSE, missing = FALSE),
    has_tm     = if_else(coalesce(tmhmm_nhelix, 0) > 0, TRUE, FALSE)
  ) |>
  distinct(input_file, .keep_all = TRUE)

num_cols <- c("comprimento","mw","gravy","pi","carga",
              "aliphatic_index","boman_index","aromaticity",
              "sasa_total","frac_hydrophobic_surface","rg",
              "plddt_mean","plddt_median","frac_plddt_lt70",
              "contact_density","iupred_mean_disorder",
              "frac_ss_helix","frac_ss_sheet","frac_ss_turn",
              "n_cys","n_ss","tmhmm_nhelix")
num_cols <- intersect(num_cols, names(df))

# ==== 2) Core tables for the manuscript ====
sum_classe <- df |>
  group_by(classe) |>
  summarise(n = n(),
            across(all_of(num_cols),
                   list(mediana=~median(.x, na.rm=TRUE), iqr=~IQR(.x,na.rm=TRUE)),
                   .names="{.col}_{.fn}")) |>
  arrange(desc(n))
write.csv(sum_classe, file.path(outroot,"tables","01_resumo_por_classe.csv"), row.names = FALSE)

sum_taxon <- df |>
  group_by(taxon_guess) |>
  summarise(n = n(),
            across(all_of(num_cols),
                   list(mediana=~median(.x, na.rm=TRUE), iqr=~IQR(.x,na.rm=TRUE)),
                   .names="{.col}_{.fn}")) |>
  arrange(desc(n))
write.csv(sum_taxon, file.path(outroot,"tables","02_resumo_por_taxon.csv"), row.names = FALSE)

# Minimal atlas table (note: use all_of() for 'pi')
cols_atlas <- c("input_file","classe","taxon_guess","comprimento","n_cys","n_ss",
                "gravy","pi","carga","aliphatic_index","boman_index","aromaticity",
                "sasa_total","frac_hydrophobic_surface","rg",
                "plddt_mean","plddt_median","frac_plddt_lt70","contact_density",
                "iupred_mean_disorder","has_signal","has_tm")
cols_atlas <- intersect(cols_atlas, names(df))
atlas_min <- df %>% dplyr::select(dplyr::all_of(cols_atlas))
write.csv(atlas_min, file.path(outroot,"tables","15_atlas_table_minimal.csv"), row.names = FALSE)

# ==== 3) Exploratory plots (single main PDF) ====
pdf(file.path(outroot,"figs","03_exploratoria.pdf"), width=8, height=6)
if ("comprimento" %in% names(df)) print(
  ggplot(df, aes(comprimento, fill = classe)) +
    geom_density(alpha=.35) + theme_bw() +
    labs(title="Lenght Distribution per class", x="aa", y="density")
)
if ("n_cys" %in% names(df)) print(
  ggplot(df, aes(n_cys, fill=classe)) +
    geom_histogram(bins=25, alpha=.55, position="identity") +
    theme_bw() + labs(title="#Cys per class", x="#Cys")
)
if ("n_ss" %in% names(df)) print(
  ggplot(df, aes(n_ss, fill=classe)) +
    geom_histogram(binwidth=1, alpha=.55, position="identity") +
    theme_bw() + labs(title="Disulfide bond per class", x="#pair S–S")
)
if ("plddt_mean" %in% names(df)) print(
  ggplot(df, aes(x=classe, y=plddt_mean, fill=classe)) +
    geom_violin(trim=FALSE) + geom_boxplot(width=.15, outlier.size=.6) +
    theme_bw() + theme(axis.text.x = element_text(angle=35, hjust=1)) +
    labs(title="Structural Quality per class (pLDDT)", x=NULL, y="pLDDT (mean)")
)
if (all(c("sasa_total","gravy") %in% names(df))) print(
  ggplot(df, aes(sasa_total, gravy, color=classe)) +
    geom_point(alpha=.6, size=1) + theme_bw() +
    labs(title="SASA vs GRAVY", x="total SASA (Å²)", y="GRAVY")
)
if (all(c("frac_ss_helix","frac_ss_sheet") %in% names(df))) print(
  ggplot(df, aes(frac_ss_helix, frac_ss_sheet, color=classe)) +
    geom_point(alpha=.6, size=1) + theme_bw() +
    labs(title="SSE fraction (helix vs sheet)", x="Helix (fraction)", y="Sheet (fraction)")
)
dev.off()

# ==== 4) Curated Spearman correlation → PDF ====
num_cols <- intersect(num_cols, names(df))
M <- df[, num_cols, drop = FALSE]
M[] <- lapply(M, function(x) suppressWarnings(as.numeric(x)))
ok_vars <- sapply(M, function(x) sum(is.finite(x)) >= 3 && sd(x, na.rm=TRUE) > 0)
M2 <- M[, ok_vars, drop=FALSE]
cor_spear <- suppressWarnings(cor(M2, method="spearman", use="pairwise.complete.obs"))
cor_spear[!is.finite(cor_spear)] <- 0; diag(cor_spear) <- 1

# Keep the most variable features for plotting (or all of them if fewer are available)
var_rank <- order(apply(M2, 2, function(x) stats::var(x, na.rm=TRUE)), decreasing = TRUE)
keep_n <- min(100, ncol(M2))
keep_vars <- colnames(M2)[var_rank[seq_len(keep_n)]]
cor_cur <- cor_spear[keep_vars, keep_vars]

ord <- hclust(as.dist(1 - cor_cur), method = "average")
ord_idx <- ord$order; cor_ord <- cor_cur[ord_idx, ord_idx]

pdf(file.path(outroot,"figs","04_cor_spearman_curado.pdf"), width=10, height=8, onefile=FALSE)
ggcorrplot(cor_ord, type="lower", lab=TRUE, lab_size=3, tl.cex=10, hc.order=FALSE,
           title="Spearman Correlation")
dev.off()

# Logs and full matrices (supplementary)
writeLines(colnames(M)[!ok_vars], file.path(outroot,"supplementary","04a_vars_removidas_por_qc.txt"))
write.csv(cor_spear, file.path(outroot,"supplementary","04b_cor_spearman_matriz_completa.csv"), row.names = TRUE)

# ==== 5) PCA with ellipses → PDF ====
cols_ok <- intersect(num_cols, names(df))
X <- df[, cols_ok, drop = FALSE]
X[] <- lapply(X, function(x) suppressWarnings(as.numeric(x)))
is_ok_col <- sapply(X, function(x){ x<-x[is.finite(x)]; length(x)>=3 && length(unique(x))>1 && sd(x,na.rm=TRUE)>0 })
X2 <- X[, is_ok_col, drop=FALSE] |> tidyr::drop_na()
stopifnot(ncol(X2) >= 2, nrow(X2) >= 3)

pca <- prcomp(X2, center = TRUE, scale. = TRUE)
cols_meta <- c("input_file", intersect(c("classe","taxon_guess","familia","taxon"), names(df)))
meta_align <- df[ stats::complete.cases(X2), cols_meta, drop = FALSE ]
scores <- as_tibble(pca$x) %>% bind_cols(meta_align)

col_grp <- if ("classe" %in% names(scores)) "classe" else if ("familia" %in% names(scores)) "familia" else if ("taxon_guess" %in% names(scores)) "taxon_guess" else "input_file"
scores[[col_grp]] <- as.factor(scores[[col_grp]])

pdf(file.path(outroot,"figs","08_pca_scree.pdf"), width=10, height=7, onefile=FALSE)
factoextra::fviz_eig(pca, addlabels=TRUE)
dev.off()

pdf(file.path(outroot,"figs","09_pca_biplot.pdf"), width=12, height=9, onefile=FALSE)
factoextra::fviz_pca_biplot(pca, repel=TRUE, col.var="firebrick", col.ind="gray30")
dev.off()

pdf(file.path(outroot,"figs","09_pca_biplot_main.pdf"), width=12, height=9, onefile=FALSE)
factoextra::fviz_pca_biplot(pca, repel=TRUE, col.var="firebrick", col.ind="white")
dev.off()

pdf(file.path(outroot,"figs","10_pca_pc1_pc2_por_grupo_elipses.pdf"), width=12, height=9, onefile=FALSE)
print(
  ggplot(scores, aes(PC1, PC2, color=.data[[col_grp]], fill=.data[[col_grp]])) +
    geom_point(alpha=.5, size=1) +
    stat_ellipse(aes(group=.data[[col_grp]]), type="t", level=0.68,
                 alpha=0.18, geom="polygon", linewidth=0.4) +
    theme_bw() +
    labs(title="PCA – PC1 vs PC2 per group", color="Grupo", fill="Grupo")
)
dev.off()

# ==== 6) UMAP + HDBSCAN (single PDF) ====
X_umap <- X2; N <- nrow(X_umap)
if (N >= 3) {
  k <- max(2L, min(15L, N-1L))
  has_dups <- any(duplicated(as.data.frame(round(X_umap, 10))))
  if (has_dups) X_umap <- as.matrix(X_umap) + matrix(rnorm(N*ncol(X_umap), sd=1e-8), nrow=N)
  
  umap_xy <- uwot::umap(X_umap, n_neighbors=k, min_dist=0.1, metric="cosine", n_threads=0)
  
  umap_df <- tibble::as_tibble(umap_xy, .name_repair=~c("UMAP1","UMAP2")) %>%
    dplyr::bind_cols(
      df %>%
        dplyr::filter(stats::complete.cases(X2)) %>%
        dplyr::select(input_file, tidyselect::any_of(c("classe","familia","taxon_guess")))
    )
  
  min_pts <- max(5L, round(sqrt(N)))
  hdb <- dbscan::hdbscan(umap_xy, minPts=min_pts)
  umap_df$cluster <- factor(ifelse(hdb$cluster==0, NA, hdb$cluster))
  
  col_grp_umap <- if ("classe" %in% names(umap_df)) "classe" else if ("familia" %in% names(umap_df)) "familia" else if ("taxon_guess" %in% names(umap_df)) "taxon_guess" else "input_file"
  
  pdf(file.path(outroot,"figs","08_09_umap_figs.pdf"), width=12, height=9, onefile=TRUE)
  print( ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color=cluster)) +
           ggplot2::geom_point(alpha=.8, size=1) + ggplot2::theme_bw() +
           ggplot2::labs(title=sprintf("UMAP + HDBSCAN (k=%d, minPts=%d)", k, min_pts)) )
  print( ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color=.data[[col_grp_umap]])) +
           ggplot2::geom_point(alpha=.8, size=1) + ggplot2::theme_bw() +
           ggplot2::labs(title="UMAP per group", color="Grupo") )
  dev.off()
  
  writeLines(sprintf("N=%d; n_neighbors=%d; minPts=%d; dups=%s", N, k, min_pts, has_dups),
             file.path(outroot,"supplementary","umap_hdbscan_params.txt"))
}

# ==== 7) PERMANOVA (supplementary) ====
keep_idx <- complete.cases(X2)
X_perm <- X2[keep_idx,,drop=FALSE]; D <- dist(scale(X_perm), method="euclidean")
cols_meta <- intersect(c("classe","taxon_guess"), names(df))
meta_perm <- df[ stats::complete.cases(X2), cols_meta, drop = FALSE ]
run_permanova <- function(D, meta, var, out_path){
  if (!(var %in% names(meta))) return(invisible(NULL))
  g <- factor(meta[[var]]); tb <- table(g, useNA="ifany")
  keep_levels <- names(tb)[tb >= 2]
  if (length(keep_levels) < 2) {
    writeLines(c(sprintf("PERMANOVA (%s): <2 levels with n>=2", var), capture.output(print(tb))),
               file.path(out_path, sprintf("permanova_%s_NOT_RUN.txt", var)))
    return(invisible(NULL))
  }
  idx <- g %in% keep_levels
  fit <- vegan::adonis2(as.dist(as.matrix(D)[idx,idx]) ~ factor(g[idx]), permutations=999)
  capture.output(fit, file=file.path(out_path, sprintf("permanova_%s.txt", var)))
  write.csv(as.data.frame(tb), file.path(out_path, sprintf("permanova_%s_counts.csv", var)))
}
run_permanova(D, meta_perm, "classe", file.path(outroot,"supplementary"))
run_permanova(D, meta_perm, "taxon_guess", file.path(outroot,"supplementary"))

# ==== 8) Between-group tests (supplementary) ====
if ("plddt_mean" %in% names(df)) {
  kr_plddt <- df |> filter(!is.na(plddt_mean)) |> kruskal_test(plddt_mean ~ classe)
  write.csv(kr_plddt, file.path(outroot,"supplementary","12_kruskal_plddt_por_classe.csv"), row.names=FALSE)
}
if ("sasa_total" %in% names(df)) {
  kr_sasa <- df |> filter(!is.na(sasa_total)) |> kruskal_test(sasa_total ~ classe)
  write.csv(kr_sasa, file.path(outroot,"supplementary","13_kruskal_sasa_por_classe.csv"), row.names=FALSE)
}

# ==== 9) Optional compact heatmap (PDF) ====
mat_scaled <- scale(X2)  # X2: amostras x variáveis (já limpo)
n_draw <- min(300, nrow(mat_scaled))
set.seed(1)
idx <- sample(seq_len(nrow(mat_scaled)), size = n_draw)

mat_sub <- mat_scaled[idx, , drop = FALSE]

# Sample-level row annotation aligned to the same subset
ann <- df %>%
  dplyr::filter(stats::complete.cases(X2)) %>%
  dplyr::slice(idx) %>%
  dplyr::transmute(Classe = factor(classe),
                   Taxon  = factor(taxon_guess)) %>%
  as.data.frame()

# Keep row names consistent between the matrix and annotation data frame
rownames(mat_sub) <- paste0("s", idx)
rownames(ann)     <- rownames(mat_sub)

# Row annotation (not top_annotation)
row_ann <- ComplexHeatmap::rowAnnotation(df = ann)

pdf(file.path(outroot, "figs", "14_heatmap_features_curto.pdf"),
    width = 12, height = 9, onefile = FALSE)
ht <- ComplexHeatmap::Heatmap(mat_sub, name = "z",
                              show_row_names = FALSE,
                              left_annotation = row_ann)
ComplexHeatmap::draw(ht)
dev.off()

# ==== 10) Optional clustering on PCs (supplementary) ====
pcs_use <- paste0("PC", 1:min(10, ncol(scores)-0))
Z <- as.matrix(scores[, pcs_use, drop=FALSE])

# k-means baseline (k = 6)
set.seed(123); k <- 6
km <- kmeans(Z, centers=k, nstart=50); scores$cluster_km <- factor(km$cluster)
pdf(file.path(outroot,"figs","11_pca_kmeans.pdf"), width=12, height=9)
print( ggplot(scores, aes(PC1, PC2, color=cluster_km)) + geom_point(alpha=.7, size=1) +
         theme_bw() + labs(title=sprintf("PCA – k-means (k=%d)", k), color="k-means") )
dev.off()

# HDBSCAN with minPts scaled to sample size
min_pts <- max(5L, round(sqrt(nrow(Z))))
hdb <- dbscan::hdbscan(Z, minPts=min_pts); scores$cluster_hdb <- factor(ifelse(hdb$cluster==0, NA, hdb$cluster))
pdf(file.path(outroot,"figs","12_pca_hdbscan.pdf"), width=12, height=9)
print( ggplot(scores, aes(PC1, PC2, color=cluster_hdb)) + geom_point(alpha=.7, size=1) +
         theme_bw() + labs(title=sprintf("PCA – HDBSCAN (minPts=%d)", min_pts), color="HDBSCAN") )
dev.off()

# Supervised LDA, if at least two levels are present
if (length(levels(scores[[col_grp]])) >= 2) {
  lda_fit <- MASS::lda(x=Z, grouping=scores[[col_grp]])
  lda_sc  <- as.data.frame(predict(lda_fit)$x); lda_sc[[col_grp]] <- scores[[col_grp]]
  pdf(file.path(outroot,"figs","15_lda_pcs.pdf"), width=12, height=9)
  print( ggplot(lda_sc, aes(LD1, LD2, color=.data[[col_grp]])) + geom_point(alpha=.8, size=1.2) +
           theme_bw() + labs(title="LDA for PCs", color="Grupo") )
  dev.off()
}

# Random forest PC importance (supplementary table)
set.seed(123)
rf <- randomForest(x=Z, y=scores[[col_grp]], ntree=500)
write.csv(randomForest::importance(rf), file.path(outroot,"supplementary","16_rf_importancia_pcs.csv"))


####10####
infile <- "attrubutes_prot.txt"
outdir <- "atlas_out2"
# -------------------------
# Packages
# -------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(ggplot2)
  library(vegan)      # PERMANOVA (adonis2)
  library(MASS)       # LDA
  library(ranger)     # Random Forest rápido
  library(scales)     # percent()
})

# -------------------------
# Helper functions
# -------------------------
write_tbl <- function(x, file) {
  readr::write_csv(x, file)
  message("  - saved: ", file)
}

find_col <- function(nms, patterns, prefer = NULL) {
  nms_low <- tolower(nms)
  hits <- unique(unlist(lapply(patterns, function(p) which(grepl(p, nms_low)))))
  if (!is.null(prefer)) {
    prefer_low <- tolower(prefer)
    hits <- c(which(nms_low %in% prefer_low), setdiff(hits, which(nms_low %in% prefer_low)))
  }
  if (length(hits) > 0) nms[hits[1]] else NA_character_
}

# -------------------------
# Input (RStudio-friendly)
# -------------------------
if (is.na(infile)) {
  message("Select the attribute file (TSV/CSV/TXT):")
  infile <- file.choose()
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Try to detect the input delimiter
sep_guess <- if (grepl("\\.csv$", infile, ignore.case = TRUE)) "," else "\t"
df_raw <- tryCatch({
  suppressMessages(readr::read_delim(infile, delim = sep_guess, guess_max = 1e5))
}, error = function(e) {
  suppressMessages(readr::read_table(infile, guess_max = 1e5))
})
df <- janitor::clean_names(df_raw)
stopifnot(nrow(df) > 0)

# -------------------------
# Sample ID and grouping
# -------------------------
id_col <- find_col(names(df), c("^id$", "input[_-]?file", "seq.*id", "name", "header"))
if (is.na(id_col)) {
  df <- df %>% mutate(sample_id = paste0("row", dplyr::row_number()))
  id_col <- "sample_id"
}
group_col <- find_col(names(df), c("^group$", "classe", "class", "label", "family", "taxon"))
if (is.na(group_col)) {
  # Use the base filename, remove the extension, and keep the prefix before the first '_'
  base_ids <- basename(as.character(df[[id_col]]))
  base_ids <- sub("\\..*$", "", base_ids)          # remove file extension
  groups   <- sub("_.*$", "", base_ids)            # keep the prefix before the first '_'
  df$group <- groups
  group_col <- "group"
} else {
  names(df)[names(df) == group_col] <- "group"
  group_col <- "group"
}

# -------------------------
# Identify key columns
# -------------------------
col_plddt_mean    <- find_col(names(df), c("plddt.*mean", "^mean.*plddt", "^plddt$"))
col_plddt_median  <- find_col(names(df), c("plddt.*median", "median.*plddt"))
col_frac_lt70     <- find_col(names(df), c("frac.*plddt.*lt.*70", "plddt.*lt.*70", "frac_lt70", "lt70"))

col_signalp       <- find_col(names(df), c("^signalp$", "signal[_-]?p", "has.*signalp"))
col_signalp_clev  <- find_col(names(df), c("signalp.*cleav", "cleav.*site"))
col_tmhmm_nh      <- find_col(names(df), c("tmhmm.*n.*helix", "n_?tm", "num_?tm", "tm_count"))
col_iupred_mean   <- find_col(names(df), c("iupred.*mean", "mean.*iupred", "^iupred$","disorder.*mean"))

# -------------------------
# pLDDT quality control
# -------------------------
qc_df <- df %>% mutate(
  plddt_mean = if (!is.na(col_plddt_mean)) suppressWarnings(as.numeric(.data[[col_plddt_mean]])) else NA_real_,
  plddt_median = if (!is.na(col_plddt_median)) suppressWarnings(as.numeric(.data[[col_plddt_median]])) else NA_real_,
  frac_plddt_lt70 = if (!is.na(col_frac_lt70)) suppressWarnings(as.numeric(.data[[col_frac_lt70]])) else NA_real_
) %>%
  mutate(
    qc_3d_exclude = ifelse(!is.na(frac_plddt_lt70) & frac_plddt_lt70 > 0.5, TRUE, FALSE),
    qc_flag = dplyr::case_when(
      qc_3d_exclude ~ "EXCLUDE_3D",
      !is.na(plddt_mean) & plddt_mean < 70 ~ "LOW_CONF_WARN",
      TRUE ~ "OK"
    )
  )

qc_summary <- qc_df %>% count(qc_flag, name = "n") %>% mutate(prop = n / sum(n))
write_tbl(qc_summary, file.path(outdir, "qc_summary.csv"))

pdf(file.path(outdir, "qc_plddt_hist.pdf"), width = 7, height = 5)
if (!all(is.na(qc_df$plddt_mean))) {
  print(
    ggplot(qc_df, aes(x = plddt_mean)) +
      geom_histogram(bins = 30) +
      labs(title = "pLDDT distribution (mean)", x = "pLDDT mean", y = "N") +
      theme_minimal()
  )
} else {
  plot.new(); title("pLDDT not available in the input file")
}
dev.off()

# -------------------------
# SignalP, TMHMM, and IUPred summaries
# -------------------------
percent_list <- list()

if (!is.na(col_signalp)) {
  sig_raw <- as.character(df[[col_signalp]])
  sig_std <- dplyr::case_when(
    tolower(sig_raw) %in% c("y","yes","true","1","signalpeptide","sp","positive","+") ~ "Positive",
    tolower(sig_raw) %in% c("n","no","false","0","none","negative","-") ~ "Negative",
    TRUE ~ ifelse(is.na(sig_raw), "NA", "Other")
  )
  percent_list$SignalP <- tibble::tibble(
    metric = "SignalP",
    level = c("Positive","Negative","Other/NA"),
    n = c(sum(sig_std=="Positive"), sum(sig_std=="Negative"), sum(sig_std %in% c("Other","NA")))
  ) %>% mutate(perc = n/sum(n))
}

if (!is.na(col_signalp_clev)) {
  clev_raw <- as.character(df[[col_signalp_clev]])
  clev_std <- ifelse(is.na(clev_raw) | clev_raw=="", "NoSite/NA", "HasCleavage")
  percent_list$SignalP_cleavage <- tibble::tibble(
    metric="SignalP_cleavage",
    level=c("HasCleavage","NoSite/NA"),
    n=c(sum(clev_std=="HasCleavage"), sum(clev_std=="NoSite/NA"))
  ) %>% mutate(perc = n/sum(n))
}

if (!is.na(col_tmhmm_nh)) {
  tm_n <- suppressWarnings(as.numeric(df[[col_tmhmm_nh]]))
  percent_list$TMHMM <- tibble::tibble(
    metric="TMHMM_n_helix",
    level=c("0 TM",">=1 TM","NA"),
    n=c(sum(tm_n==0, na.rm=TRUE), sum(tm_n>=1, na.rm=TRUE), sum(is.na(tm_n)))
  ) %>% mutate(perc = n/sum(n))
}

if (!is.na(col_iupred_mean)) {
  dis <- suppressWarnings(as.numeric(df[[col_iupred_mean]]))
  cat_iup <- dplyr::case_when(
    is.na(dis) ~ "NA",
    dis <= 0.2 ~ "Structured (≤0.2)",
    dis <= 0.5 ~ "Intermediate (0.2–0.5)",
    TRUE ~ "Disordered (>0.5)"
  )
  percent_list$IUPred <- tibble::tibble(metric="IUPred_mean", level=names(sort(table(cat_iup), TRUE))) %>%
    dplyr::left_join(tibble::tibble(level=unique(cat_iup), n=as.numeric(table(cat_iup))), by="level") %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n), perc = n/sum(n))
}

if (length(percent_list) > 0) {
  percent_df <- dplyr::bind_rows(percent_list)
  write_tbl(percent_df, file.path(outdir, "signalp_tmhmm_iupred_summary.csv"))
} else {
  message("No SignalP/TMHMM/IUPred columns detected for summary.")
}

# -------------------------
# PCA after removing EXCLUDE_3D entries
# -------------------------
df_use <- qc_df %>% dplyr::filter(!qc_3d_exclude)

num_cols <- names(df_use)[sapply(df_use, is.numeric)]
num_cols <- setdiff(num_cols, c("qc_3d_exclude"))
stopifnot(length(num_cols) >= 2)

X <- df_use %>% dplyr::select(dplyr::all_of(num_cols)) %>% as.data.frame()
X_scaled <- scale(X)

pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
var_tbl <- tibble::tibble(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  cum_variance = cumsum(var_exp)
)
write_tbl(var_tbl, file.path(outdir, "pca_variance.csv"))

scores <- as.data.frame(pca$x) %>%
  dplyr::mutate(!!id_col := df_use[[id_col]], group = df_use[[group_col]])
write_tbl(scores, file.path(outdir, "pca_scores.csv"))

pc1ve <- scales::percent(var_exp[1]); pc2ve <- scales::percent(var_exp[2])
pdf(file.path(outdir, "pca_scatter.pdf"), width = 7, height = 5)
print(
  ggplot(scores, aes(x = PC1, y = PC2, color = factor(group))) +
    geom_point(size = 2, alpha = 0.85) +
    labs(title = "PCA",
         x = paste0("PC1 (", pc1ve, ")"),
         y = paste0("PC2 (", pc2ve, ")"),
         color = "Group") +
    theme_minimal()
)
dev.off()

# -------------------------
# PERMANOVA
# -------------------------
permanova_msg <- "PERMANOVA not run (insufficient groups or too few observations)."
try({
  if (!is.null(scores$group) && length(unique(scores$group)) >= 2) {
    grp <- factor(scores$group)
    X_scaled_use <- as.data.frame(X_scaled)[match(df_use[[id_col]], df_use[[id_col]]), , drop=FALSE]
    ok_rows <- stats::complete.cases(X_scaled_use) & !is.na(grp)
    X_scaled_ok <- X_scaled_use[ok_rows, , drop=FALSE]
    grp_ok <- droplevels(grp[ok_rows])
    
    if (nrow(X_scaled_ok) >= 10 && length(unique(grp_ok)) >= 2) {
      set.seed(123)
      dist_mat <- dist(X_scaled_ok, method = "euclidean")
      perma <- vegan::adonis2(dist_mat ~ grp_ok, permutations = 999)
      capture.output(perma, file = file.path(outdir, "permanova_result.txt"))
      permanova_msg <- NULL
    }
  }
}, silent = TRUE)
if (!is.null(permanova_msg)) writeLines(permanova_msg, con = file.path(outdir, "permanova_result.txt"))

# -------------------------
# Lightweight classifiers on PC scores
# -------------------------
k_auto <- max(2, which(cumsum(var_exp) <= 0.80))
k_auto <- min(max(k_auto, 2), min(ncol(scores) - 2, 20))
pc_cols <- paste0("PC", seq_len(k_auto))

clf_df <- scores %>%
  dplyr::select(dplyr::all_of(c(pc_cols, id_col, "group"))) %>%
  tidyr::drop_na(dplyr::all_of(c(pc_cols, "group")))

if (nrow(clf_df) >= 20 && length(unique(clf_df$group)) >= 2) {
  set.seed(42)
  idx <- unlist(tapply(seq_len(nrow(clf_df)), clf_df$group, function(ix) {
    n <- length(ix); sample(ix, size = max(1, floor(0.8 * n)))
  }))
  train <- clf_df[idx, , drop=FALSE]
  test  <- clf_df[-idx, , drop=FALSE]
  
  acc_res <- list()
  
  lda_fit <- tryCatch(MASS::lda(group ~ ., data = dplyr::select(train, -dplyr::all_of(id_col))), error = function(e) NULL)
  if (!is.null(lda_fit)) {
    lda_pred <- predict(lda_fit, newdata = dplyr::select(test, -dplyr::all_of(c(id_col, "group"))))$class
    acc_res$LDA <- tibble::tibble(model="LDA", k_pcs = k_auto, accuracy_holdout = mean(lda_pred == test$group))
  }
  
  rf_fit <- tryCatch(
    ranger::ranger(group ~ ., data = dplyr::select(train, -dplyr::all_of(id_col)),
                   probability = FALSE, num.trees = 500, mtry = floor(sqrt(k_auto)),
                   min.node.size = 3, respect.unordered.factors = TRUE, seed = 42),
    error = function(e) NULL
  )
  if (!is.null(rf_fit)) {
    rf_pred <- predict(rf_fit, data = dplyr::select(test, -dplyr::all_of(c(id_col, "group"))))$predictions
    acc_res$RandomForest <- tibble::tibble(model="RandomForest", k_pcs = k_auto, accuracy_holdout = mean(rf_pred == test$group))
  }
  
  if (length(acc_res) > 0) {
    acc_tab <- dplyr::bind_rows(acc_res)
    write_tbl(acc_tab, file.path(outdir, "classifier_accuracy.csv"))
  } else {
    writeLines("Classifiers could not be trained (collinearity or a single group).",
               con = file.path(outdir, "classifier_accuracy.txt"))
  }
} else {
  writeLines("Classifiers not run (insufficient samples or groups).",
             con = file.path(outdir, "classifier_accuracy.txt"))
}

# -------------------------
# Note on per-helix hydrophobic moment (DSSP)
# -------------------------
note_txt <- c(
  "Note: per-helix hydrophobic moment (μH, DSSP-based) was not calculated here.",
  "This requires residue-level DSSP output (or helix boundaries) for each sequence.",
  "Once residue-level secondary structure data are available (or PDB/mmCIF files to generate DSSP), ",
  "the per-helix μH calculation can be added and compared with SASA and charge."
)
writeLines(note_txt, con = file.path(outdir, "note_muH_needed.txt"))

message("\nDone. Outputs written to: ", normalizePath(outdir))

gc()
