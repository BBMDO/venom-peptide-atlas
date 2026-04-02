#A Structure-Informed Atlas of Venom-Derived Peptides Reveals the Organization of Chemical Space
#Thaís Caroline Gonçalves, Eduardo Henrique Toral Cortez, and Danilo Trabuco Amaral

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(uwot)
  library(randomForest)
  library(ggrepel)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)

venom_file   <- ifelse(length(args) >= 1, args[1], "attrubutes_prot_venom.tsv")
nonven_file  <- ifelse(length(args) >= 2, args[2], "param_non-venom.tsv")
outdir       <- ifelse(length(args) >= 3, args[3], "results_venom_vs_nonvenom")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# =========================
# 1. Flexible input reader
# =========================
read_table_flexible <- function(file) {
  first_line <- readLines(file, n = 1, warn = FALSE)
  
  tab_count   <- stringr::str_count(first_line, "\t")
  comma_count <- stringr::str_count(first_line, ",")
  semi_count  <- stringr::str_count(first_line, ";")
  
  if (tab_count >= comma_count && tab_count >= semi_count) {
    df <- readr::read_tsv(file, show_col_types = FALSE, progress = FALSE)
  } else if (comma_count >= semi_count) {
    df <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
  } else {
    df <- readr::read_delim(file, delim = ";", show_col_types = FALSE, progress = FALSE)
  }
  
  as_tibble(df)
}

# =========================
# 2. Basic column name standardization
# =========================
clean_names_simple <- function(df) {
  colnames(df) <- colnames(df) %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[[:space:]]+", "_") %>%
    stringr::str_replace_all("[^[:alnum:]_]", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
  df
}

venom_raw  <- read_table_flexible(venom_file)  %>% clean_names_simple()
nonven_raw <- read_table_flexible(nonven_file) %>% clean_names_simple()

cat("Venom dimensions:", dim(venom_raw), "\n")
cat("Non-venom dimensions:", dim(nonven_raw), "\n")

# =========================
# 3. Shared feature space
# =========================
common_cols <- intersect(colnames(venom_raw), colnames(nonven_raw))
common_cols <- setdiff(common_cols, c("dataset", "class"))

if (length(common_cols) < 5) {
  stop("Too few shared columns were detected after column-name standardization.")
}

writeLines(common_cols, file.path(outdir, "common_columns.txt"))

venom_shared <- venom_raw %>%
  dplyr::select(all_of(common_cols)) %>%
  mutate(dataset = "venom")

nonven_shared <- nonven_raw %>%
  dplyr::select(all_of(common_cols)) %>%
  mutate(dataset = "non_venom")

all_shared <- bind_rows(venom_shared, nonven_shared)

write_tsv(all_shared, file.path(outdir, "all_shared_space.tsv"))

# =========================
# 4. Safe numeric conversion
# =========================
safe_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x2 <- suppressWarnings(as.numeric(as.character(x)))
  x2
}

candidate_numeric <- setdiff(colnames(all_shared), "dataset")

all_shared_num <- all_shared %>%
  mutate(across(all_of(candidate_numeric), safe_as_numeric))

numeric_cols <- all_shared_num %>%
  dplyr::select(-dataset) %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# Remove columns containing only missing values
numeric_cols <- numeric_cols[colSums(!is.na(all_shared_num[, numeric_cols, drop = FALSE])) > 0]

# Remove zero-variance columns
numeric_cols <- numeric_cols[sapply(all_shared_num[, numeric_cols, drop = FALSE], function(x) {
  x <- x[!is.na(x)]
  length(unique(x)) > 1
})]

if (length(numeric_cols) < 3) {
  stop("Fewer than three informative shared numeric columns are available for downstream analyses.")
}

# =========================
# 5. Random equal-N balancing
# =========================
venom_df  <- all_shared_num %>% filter(dataset == "venom")
nonven_df <- all_shared_num %>% filter(dataset == "non_venom")

venom_n  <- nrow(venom_df)
nonven_n <- nrow(nonven_df)

n_sample <- min(venom_n, nonven_n)

venom_random_equalN <- venom_df %>%
  slice_sample(n = n_sample)

nonven_random_equalN <- nonven_df %>%
  slice_sample(n = n_sample)

balanced_random <- bind_rows(venom_random_equalN, nonven_random_equalN)

write_tsv(balanced_random, file.path(outdir, "balanced_random_equalN.tsv"))

# =========================
# 6. Length-matched balancing
# =========================
balanced_length <- NULL

if ("length" %in% colnames(all_shared_num)) {
  venom_len  <- venom_df %>% filter(!is.na(length))
  nonven_len <- nonven_df %>% filter(!is.na(length))
  
  if (nrow(venom_len) > 20 && nrow(nonven_len) > 20) {
    venom_len <- venom_len %>%
      mutate(row_id = row_number())
    
    nonven_len <- nonven_len %>%
      mutate(row_id = row_number())
    
    matched_list <- lapply(seq_len(nrow(venom_len)), function(i) {
      vlen <- venom_len$length[i]
      
      candidates <- nonven_len %>%
        mutate(diff_len = abs(length - vlen)) %>%
        arrange(diff_len)
      
      if (nrow(candidates) == 0) return(NULL)
      
      best <- candidates[1, ]
      nonven_len <<- nonven_len %>% filter(row_id != best$row_id)
      
      tibble(
        venom_row_id = venom_len$row_id[i],
        nonven_row_id = best$row_id,
        venom_length = vlen,
        nonven_length = best$length
      )
    })
    
    matched_pairs <- bind_rows(matched_list)
    
    if (nrow(matched_pairs) > 0) {
      venom_matched <- venom_len %>%
        filter(row_id %in% matched_pairs$venom_row_id) %>%
        select(-row_id)
      
      nonven_matched <- nonven_raw %>% 
        clean_names_simple() %>%
        dplyr::select(all_of(common_cols)) %>%
        mutate(dataset = "non_venom") %>%
        mutate(across(all_of(candidate_numeric), safe_as_numeric)) %>%
        mutate(row_id = row_number()) %>%
        filter(row_id %in% matched_pairs$nonven_row_id) %>%
        select(-row_id)
      
      balanced_length <- bind_rows(venom_matched, nonven_matched)
      
      write_tsv(balanced_length, file.path(outdir, "balanced_length_matched.tsv"))
      write_tsv(matched_pairs, file.path(outdir, "length_matched_pairs.tsv"))
    }
  }
}

# =========================
# 7. Helper functions
# =========================
impute_median <- function(df, cols) {
  df2 <- df
  for (cc in cols) {
    med <- median(df2[[cc]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    df2[[cc]][is.na(df2[[cc]])] <- med
  }
  df2
}

run_univariate_tests <- function(df, cols, prefix) {
  res <- lapply(cols, function(v) {
    sub <- df %>% select(dataset, all_of(v)) %>% filter(!is.na(.data[[v]]))
    
    if (nrow(sub) < 5) return(NULL)
    if (length(unique(sub$dataset)) < 2) return(NULL)
    if (length(unique(sub[[v]])) < 2) return(NULL)
    
    wt <- tryCatch(
      wilcox.test(sub[[v]] ~ sub$dataset),
      error = function(e) NULL
    )
    
    medians <- sub %>%
      group_by(dataset) %>%
      summarise(median = median(.data[[v]], na.rm = TRUE), .groups = "drop")
    
    tibble(
      variable = v,
      p_value = ifelse(is.null(wt), NA, wt$p.value),
      median_venom = medians$median[medians$dataset == "venom"] %>% .[1],
      median_nonvenom = medians$median[medians$dataset == "non_venom"] %>% .[1]
    )
  }) %>% bind_rows()
  
  if (nrow(res) > 0) {
    res <- res %>%
      mutate(p_adj_fdr = p.adjust(p_value, method = "fdr")) %>%
      arrange(p_adj_fdr, p_value)
    
    write_tsv(res, file.path(outdir, paste0(prefix, "_univariate_stats.tsv")))
  }
  
  invisible(res)
}

run_pca <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  x <- scale(df_imp[, cols, drop = FALSE])
  pca <- prcomp(x, center = FALSE, scale. = FALSE)
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$dataset <- df_imp$dataset
  
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  p <- ggplot(scores, aes(PC1, PC2, color = dataset)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("PCA -", prefix),
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )
  
  ggsave(file.path(outdir, paste0(prefix, "_pca.pdf")), p, width = 7, height = 5)
  
  loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE]) %>%
    rownames_to_column("variable")
  write_tsv(loadings, file.path(outdir, paste0(prefix, "_pca_loadings.tsv")))
  
  invisible(list(pca = pca, scores = scores))
}

run_umap_plot <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  x <- scale(df_imp[, cols, drop = FALSE])
  
  um <- uwot::umap(x, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", verbose = FALSE)
  umdf <- as.data.frame(um)
  colnames(umdf) <- c("UMAP1", "UMAP2")
  umdf$dataset <- df_imp$dataset
  
  p <- ggplot(umdf, aes(UMAP1, UMAP2, color = dataset)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_classic(base_size = 12) +
    labs(title = paste("UMAP -", prefix))
  
  ggsave(file.path(outdir, paste0(prefix, "_umap.pdf")), p, width = 7, height = 5)
  
  write_tsv(umdf, file.path(outdir, paste0(prefix, "_umap_coordinates.tsv")))
  invisible(umdf)
}

run_permanova <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  x <- scale(df_imp[, cols, drop = FALSE])
  
  ad <- vegan::adonis2(x ~ dataset, data = df_imp, method = "euclidean", permutations = 999)
  
  capture.output(ad, file = file.path(outdir, paste0(prefix, "_permanova.txt")))
  invisible(ad)
}

run_rf <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  df_imp$dataset <- as.factor(df_imp$dataset)
  
  rf <- randomForest(
    x = df_imp[, cols, drop = FALSE],
    y = df_imp$dataset,
    importance = TRUE,
    ntree = 500
  )
  
  imp <- importance(rf, type = 2) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    arrange(desc(MeanDecreaseGini))
  
  write_tsv(imp, file.path(outdir, paste0(prefix, "_rf_importance.tsv")))
  
  p <- ggplot(head(imp, 20), aes(x = reorder(variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(title = paste("Random Forest importance -", prefix), x = "", y = "MeanDecreaseGini")
  
  ggsave(file.path(outdir, paste0(prefix, "_rf_top20.pdf")), p, width = 7, height = 6)
  
  invisible(rf)
}

run_density_plots <- function(df, cols, prefix, top_n = 12) {
  vars_to_plot <- cols[1:min(top_n, length(cols))]
  
  pdf(file.path(outdir, paste0(prefix, "_density_panels.pdf")), width = 10, height = 8)
  for (v in vars_to_plot) {
    p <- ggplot(df, aes(x = .data[[v]], fill = dataset)) +
      geom_density(alpha = 0.4) +
      theme_classic(base_size = 11) +
      labs(title = v, x = v, y = "Density")
    print(p)
  }
  dev.off()
}

run_pca <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  
  # Keep numeric columns only
  xdf <- df_imp[, cols, drop = FALSE]
  xdf <- xdf[, sapply(xdf, is.numeric), drop = FALSE]
  
  # Remove zero-variance columns OU NA após imputação
  keep <- sapply(xdf, function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    length(x) > 1 && sd(x) > 0
  })
  
  xdf <- xdf[, keep, drop = FALSE]
  
  if (ncol(xdf) < 2) {
    stop("Fewer than two informative columns remain for PCA after removing constant features.")
  }
  
  x <- scale(xdf, center = TRUE, scale = TRUE)
  
  # Remove columns that remain problematic after scaling
  x <- x[, colSums(is.na(x)) == 0, drop = FALSE]
  
  if (ncol(x) < 2) {
    stop("Fewer than two valid columns remain for PCA after scaling.")
  }
  
  pca <- prcomp(x, center = FALSE, scale. = FALSE)
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$dataset <- df_imp$dataset
  
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  p <- ggplot(scores, aes(PC1, PC2, color = dataset)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("PCA -", prefix),
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )
  
  ggsave(file.path(outdir, paste0(prefix, "_pca.pdf")), p, width = 7, height = 5)
  
  loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("variable")
  
  write_tsv(loadings, file.path(outdir, paste0(prefix, "_pca_loadings.tsv")))
  
  invisible(list(pca = pca, scores = scores, kept_columns = colnames(xdf)))
}

run_umap_plot <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  
  xdf <- df_imp[, cols, drop = FALSE]
  xdf <- xdf[, sapply(xdf, is.numeric), drop = FALSE]
  
  # Remove constant or otherwise invalid columns
  keep <- sapply(xdf, function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    length(x) > 1 && sd(x) > 0
  })
  
  xdf <- xdf[, keep, drop = FALSE]
  
  if (ncol(xdf) < 2) {
    stop("Fewer than two informative columns remain for UMAP.")
  }
  
  # Re-apply median imputation to the selected subset as a safety step
  for (cc in colnames(xdf)) {
    med <- median(xdf[[cc]], na.rm = TRUE)
    if (is.na(med) || !is.finite(med)) med <- 0
    xdf[[cc]][is.na(xdf[[cc]]) | !is.finite(xdf[[cc]])] <- med
  }
  
  x <- scale(xdf, center = TRUE, scale = TRUE)
  x <- as.data.frame(x)
  
  # Remove columns that became invalid after scaling
  keep2 <- sapply(x, function(z) all(is.finite(z)) && sd(z) > 0)
  x <- x[, keep2, drop = FALSE]
  
  if (ncol(x) < 2) {
    stop("Fewer than two valid columns remain for UMAP after scaling.")
  }
  
  # Remove rows with any residual invalid values
  good_rows <- complete.cases(x) & apply(x, 1, function(z) all(is.finite(z)))
  x <- x[good_rows, , drop = FALSE]
  meta <- df_imp[good_rows, , drop = FALSE]
  
  if (nrow(x) < 5) {
    stop("Too few valid rows remain for UMAP.")
  }
  
  um <- uwot::umap(
    as.matrix(x),
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "euclidean",
    verbose = FALSE
  )
  
  umdf <- as.data.frame(um)
  colnames(umdf) <- c("UMAP1", "UMAP2")
  umdf$dataset <- meta$dataset
  
  p <- ggplot(umdf, aes(UMAP1, UMAP2, color = dataset)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_classic(base_size = 12) +
    labs(title = paste("UMAP -", prefix))
  
  ggsave(file.path(outdir, paste0(prefix, "_umap.pdf")), p, width = 7, height = 5)
  write_tsv(umdf, file.path(outdir, paste0(prefix, "_umap_coordinates.tsv")))
  
  invisible(umdf)
}

run_permanova <- function(df, cols, prefix) {
  df_imp <- impute_median(df, cols)
  
  xdf <- df_imp[, cols, drop = FALSE]
  xdf <- xdf[, sapply(xdf, is.numeric), drop = FALSE]
  
  # Filter problematic columns
  keep <- sapply(xdf, function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    length(x) > 1 && sd(x) > 0
  })
  
  xdf <- xdf[, keep, drop = FALSE]
  
  if (ncol(xdf) < 2) {
    stop("Fewer than two informative columns remain for PERMANOVA.")
  }
  
  # Additional imputation safeguard
  for (cc in colnames(xdf)) {
    med <- median(xdf[[cc]], na.rm = TRUE)
    if (is.na(med) || !is.finite(med)) med <- 0
    xdf[[cc]][is.na(xdf[[cc]]) | !is.finite(xdf[[cc]])] <- med
  }
  
  # Standardization
  x <- scale(xdf, center = TRUE, scale = TRUE)
  x <- as.data.frame(x)
  
  # Remove columns that became invalid after scaling
  keep2 <- sapply(x, function(z) {
    z <- as.numeric(z)
    all(is.finite(z)) && sd(z) > 0
  })
  x <- x[, keep2, drop = FALSE]
  
  if (ncol(x) < 2) {
    stop("Fewer than two valid columns remain for PERMANOVA after scaling.")
  }
  
  # Remove rows with residual invalid values
  good_rows <- complete.cases(x) & apply(x, 1, function(z) all(is.finite(z)))
  x <- x[good_rows, , drop = FALSE]
  meta <- df_imp[good_rows, , drop = FALSE]
  
  if (nrow(x) < 5) {
    stop("Too few valid rows remain for PERMANOVA.")
  }
  
  ad <- vegan::adonis2(
    as.matrix(x) ~ dataset,
    data = meta,
    method = "euclidean",
    permutations = 999
  )
  
  capture.output(ad, file = file.path(outdir, paste0(prefix, "_permanova.txt")))
  invisible(ad)
}

# =========================
# 8. Run analyses
# =========================
message("Running analyses on the random equal-N dataset...")
run_univariate_tests(balanced_random, numeric_cols, "random_equalN")
run_pca(balanced_random, numeric_cols, "random_equalN")
run_umap_plot(balanced_random, numeric_cols, "random_equalN")
run_permanova(balanced_random, numeric_cols, "random_equalN")
run_rf(balanced_random, numeric_cols, "random_equalN")
run_density_plots(balanced_random, numeric_cols, "random_equalN")

if (!is.null(balanced_length)) {
  message("Running analyses on the length-matched dataset...")
  run_univariate_tests(balanced_length, numeric_cols, "length_matched")
  run_pca(balanced_length, numeric_cols, "length_matched")
  run_umap_plot(balanced_length, numeric_cols, "length_matched")
  run_permanova(balanced_length, numeric_cols, "length_matched")
  run_rf(balanced_length, numeric_cols, "length_matched")
  run_density_plots(balanced_length, numeric_cols, "length_matched")
}

# =========================
# 9. Final run summary
# =========================
summary_txt <- c(
  paste("Venom file:", venom_file),
  paste("Non-venom file:", nonven_file),
  paste("Shared columns:", length(common_cols)),
  paste("Numeric shared columns used:", length(numeric_cols)),
  paste("Venom N:", venom_n),
  paste("Non-venom N:", nonven_n),
  paste("Balanced random equal N:", nrow(balanced_random)),
  paste("Length-matched generated:", !is.null(balanced_length))
)

writeLines(summary_txt, file.path(outdir, "run_summary.txt"))

message("Analysis completed successfully.")
