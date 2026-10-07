downstreamAnalysis <- function(hvgListType, name) {
  
  dwnstrm_results <- list()
  
  variability_data <- list()
  
  gene_means <- list()
  gene_vars <- list()
  
  for (sample in names(seurat_list)) {
    obj <- seurat_list[[sample]]
    
    genes <- hvgListType[[sample]]
    
    counts <- GetAssayData(
      obj,
      layer = "counts",
      assay = "RNA"
    )[genes, ]
    
    gene_means[[sample]] <- rowMeans(counts)
    gene_vars[[sample]] <- apply(
      counts, 1, var
    )
    
    variability_data[[sample]] <- counts
  }
  
  gene_df <- data.frame(gene_means = gene_means[["L1"]], gene_vars = gene_vars[["L1"]])
  
  scatter <- ggplot(data = gene_df) +
    geom_point(mapping = aes(x = gene_means,
                             y = gene_vars)) +
    labs(title = "Mean-Variance Relationship in Sample L1",
         x = "Mean Expression",
         y = "Variance")
  
  ggsave(file.path("Graphs", name, "mean_vs_var.png"), plot = scatter)
  
  # Structure data and conduct NB modeling, retrieve dispersion results
  
  dispersion_results <- list()
  
  for (sample in names(seurat_list)) {
    obj <- seurat_list[[sample]]
    
    for (gene in hvgListType[[sample]]) {
      counts <- GetAssayData(
        obj,
        layer = "counts",
        assay = "RNA"
      )[gene, ]
      
      df <- data.frame(
        counts = as.numeric(counts)
      )
      
      fit <- glmmTMB(
        counts ~ 1,
        ziformula = ~0,
        family = nbinom2,
        data = df
      )
      
      if (fit$sdr$pdHess) {
        dispersion_results[[sample]][gene] <- sigma(fit)
      } else {
        dispersion_results[[sample]][gene] <- NA
      }
    }
  }
  
  dwnstrm_results$dispersion_results <- dispersion_results
  
  gene_disp_df <- data.frame()
  
  for (sample in names(dispersion_results)) {
    
    values <- dispersion_results[[sample]]
    
    temp <- data.frame(
      sample = sample,
      gene = names(values),
      dispersion = as.numeric(values)
    )
    
    gene_disp_df <- rbind(gene_disp_df, temp)
    
  }
  
  gene_rank <- gene_disp_df %>%
    group_by(gene) %>%
    summarise(
      across_sample_var = var(dispersion, na.rm = TRUE),
      mean_dispersion = mean(dispersion, na.rm = TRUE),
      range = max(dispersion, na.rm = TRUE) - min(dispersion, na.rm = TRUE),
      n_samples = n_distinct(sample)
    ) %>%
    filter(n_samples >= 2) %>%
    arrange(desc(across_sample_var))
  
  top_genes <- head(gene_rank, 20)
  
  tmpgraph <- ggplot(top_genes,
                     aes(x = across_sample_var,
                         y = reorder(gene, across_sample_var))) +
    geom_col(fill = "steelblue") +
    theme_minimal() +
    labs(
      title = "Genes contributing most to variability differences across samples",
      x = "Between-sample variance in dispersion",
      y = "Gene"
    )
  
  ggsave(file.path("Graphs", name, "tempGraph.png"), plot = tmpgraph)
  
  dwnstrm_results$top_genes <- top_genes
  
  # Visualizing the NB model
  
  nb_params <- data.frame()
  
  for (sample in names(seurat_list)) {
    
    genes <- names(dispersion_results[[sample]])
    
    temp <- data.frame(
      sample = sample,
      gene = genes,
      mean_expression = gene_means[[sample]][genes],
      dispersion = dispersion_results[[sample]][genes]
    )
    
    nb_params <- rbind(nb_params, temp)
  }
  
  dwnstrm_results$nb_params <- nb_params
  
  disp <- nb_params$dispersion
  disp <- disp[is.finite(disp) & disp < 10]
  
  mu <- mean(nb_params$mean_expression, na.rm = TRUE)
  alpha <- median(disp, na.rm = TRUE)
  
  size <- 1 / alpha
  prob <- size / (size + mu)
  
  dwnstrm_results$size <- size
  dwnstrm_results$prob <- prob
  
  counts_simulation <- rnbinom(
    n = 10000,
    size = size,
    prob = prob
  )
  
  nb_data <- data.frame(count = counts_simulation)
  
  dwnstrm_results$nb_data <- nb_data
  
  bar <- ggplot(nb_data, mapping = aes(x = count)) +
    geom_bar(
      fill = "slateblue",
      color = "black"
    ) +
    theme_minimal() +
    labs(title = "Negative Binomial Distribution Based on Observed Mean and Dispersion",
         x = "Feature Expression Count",
         y = "Frequency")
  
  ggsave(file.path("Graphs", name, "barplot.png"), plot = bar)
  
  # Visualizing NB model with multiple different dispersion factors
  
  disp_values <- list()
  
  for (sample in names(dispersion_results)) {
    disp_values[[sample]] <- median(dispersion_results[[sample]][!is.na(dispersion_results[[sample]])])
  }
  
  disp_values[["Test_Large"]] <- 0.25
  
  mult_disp_data <- data.frame()
  
  for (al in disp_values) {
    size <- 1 / al
    prob <- size / (size + mu)
    
    mult_counts <- rnbinom(
      n = 10000,
      size = size,
      prob = prob
    )
    
    temp <- data.frame(
      count = mult_counts,
      dispersion = as.factor(al),
      sample = names(disp_values[disp_values == al])
    )
    
    mult_disp_data <- rbind(mult_disp_data, temp)
  }
  
  dwnstrm_results$mult_disp_data <- mult_disp_data
  
  
  multbar <- ggplot(mult_disp_data, mapping = aes(x = count, fill = sample)) +
    geom_bar(
      position = "identity",
      alpha = 0.3
    ) +
    theme_minimal() +
    facet_wrap(~sample)
  
  ggsave(file.path("Graphs", name, "mult_barplot.png"), plot = multbar)
  
  # Create data frame and visualize dispersion results
  
  results_df <- data.frame()
  
  for (sample in names(dispersion_results)) {
    values <- dispersion_results[[sample]][dispersion_results[[sample]] < 10 & !is.na(dispersion_results[[sample]])]
    
    temp <- data.frame(
      sample = sample,
      dispersion = values
    )
    
    results_df <- rbind(results_df, temp)
  }
  
  results_df$sample_type <- ifelse(
    results_df$sample == "LN",
    "Normal",
    "Tumor"
  )
  
  wilcox <- wilcox.test(
    log(dispersion) ~ sample_type,
    data = results_df
  )
  
  results_df$log_dispersion <- log(results_df$dispersion)
  
  # within group variance btw tumor subgroups may have lead to smaller size effect
  
  w_effsize <- wilcox_effsize(
    data = results_df,
    formula = log_dispersion ~ sample_type
  )
  
  dwnstrm_results$results <- results_df
  dwnstrm_results$wilcox <- wilcox
  dwnstrm_results$w_effsize <- w_effsize
  
  bxplt <- ggplot(data = results_df) +
    geom_boxplot(mapping = aes(sample, -log(dispersion))) +
    labs(title = "Gene Expression Variability Dispersion Levels by Sample",
         x = "Sample Type",
         y = "-log(Dispersion Value)") +
    theme_minimal()
  
  ggsave(file.path("Graphs", name, "boxplot.png"), plot = bxplt)
  
  vln <- ggplot(data = results_df,
                mapping = aes(sample, -log(dispersion), fill = sample)) +
    geom_violin() +
    stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.75) +
    labs(title = "Gene Expression Variability Dispersion Levels by Sample",
         x = "Sample Type",
         y = "-log(Dispersion Value)")
  
  ggsave(file.path("Graphs", name, "vlnplot.png"), plot = vln)
  
  dwnstrm_results$kruskal <- kruskal.test(dispersion ~ sample, data = results_df)
  
  effSize <- rstatix::kruskal_effsize(dispersion ~ sample, data = results_df)
  
  dwnstrm_results$effSize <- effSize
  
  uspet <- upset(fromList(hvgListType), order.by = "freq", nsets = 6)
  
  # This is showing that distributions differ slightly, but 
  # HVG identity is largely non-overlapping across samples
  
  dunn_results <- dunn_test(results_df, dispersion ~ sample, p.adjust.method = "bonferroni")
  
  dwnstrm_results$dunn_results <- dunn_results
  
  dwnstrm_results$upset <- upset
  
  return(dwnstrm_results)
}