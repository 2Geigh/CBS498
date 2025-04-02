library(pheatmap)
library(gprofiler2)
library(reshape2)

# Define core parameters
CRITICAL_VALUE <- 0.01
lfc_cutoff <- 2
base_output_dir <- "./results"
heatmap_dir <- file.path(base_output_dir, "heatmaps")

# Create output directories
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Input files organized by timepoint
timepoints <- list(
  "1h" = c(
    "./deseq2_output/timepoint_1h/SigGenes_HopAB1_vs_D36E_EV_1h.csv",
    "./deseq2_output/timepoint_1h/SigGenes_HopB1_vs_D36E_EV_1h.csv",
    "./deseq2_output/timepoint_1h/SigGenes_HopN1_vs_D36E_EV_1h.csv"
  ),
  "8h" = c(
    "./deseq2_output/timepoint_8h/SigGenes_HopAB1_vs_D36E_EV_8h.csv",
    "./deseq2_output/timepoint_8h/SigGenes_HopB1_vs_D36E_EV_8h.csv",
    "./deseq2_output/timepoint_8h/SigGenes_HopN1_vs_D36E_EV_8h.csv"
  )
)

# Universal heatmap styling parameters
heatmap_style <- list(
  color_palette = colorRampPalette(c("green", "black", "red"))(50),
  fontsize_row = 6,
  fontsize_col = 9,
  angle_col = 45,
  cluster_cols = FALSE,
  clustering_method = "ward.D2",
  border_color = NA,
  show_rownames = TRUE
)

# Process each timepoint separately
for (timepoint in names(timepoints)) {
    # Select files for current timepoint
    current_files <- timepoints[[timepoint]]

    # Aggregate DEGs for this timepoint
    deg_list <- list()
    for (file in current_files) {
        deg_data <- read.csv(file)
        colnames(deg_data)[1] <- "gene_id"
        sig_genes <- deg_data[deg_data$padj < CRITICAL_VALUE & 
                                abs(deg_data$log2FoldChange) > lfc_cutoff, ]
        deg_list[[basename(file)]] <- sig_genes$gene_id
    }

    # Create timepoint-specific expression matrix
    all_degs <- unique(unlist(deg_list))
    sample_names <- gsub("SigGenes_|.csv", "", basename(current_files))
    expression_matrix <- matrix(NA, nrow = length(all_degs), ncol = length(current_files),
                                dimnames = list(all_degs, sample_names))

    # Populate matrix with log2FoldChange values
    for (i in seq_along(current_files)) {
        deg_data <- read.csv(current_files[i])
        colnames(deg_data)[1] <- "gene_id"
        matched_rows <- match(all_degs, deg_data$gene_id)
        expression_matrix[,i] <- deg_data$log2FoldChange[matched_rows]
    }
    message("Original DEG counts per sample:")
    print(colSums(!is.na(expression_matrix)))

    # Matrix cleaning with timepoint-specific parameters
    expression_matrix_clean <- expression_matrix

    # Keep genes with ≥2 observation in this timepoint (to reduce noise)
    expression_matrix_clean <- expression_matrix_clean[rowSums(!is.na(expression_matrix)) >= 2, ]

    # Winsorize using timepoint-specific distribution
    lfc_quantiles <- quantile(expression_matrix_clean, probs = c(0.01, 0.99), na.rm = TRUE)
    cap_value <- ceiling(max(abs(lfc_quantiles)) * 10)/10
    expression_matrix_clean[expression_matrix_clean < -cap_value] <- -cap_value
    expression_matrix_clean[expression_matrix_clean > cap_value] <- cap_value
    expression_matrix_clean[is.na(expression_matrix_clean)] <- 0

    # Generate heatmap
    if(nrow(expression_matrix_clean) > 1) {
        # Generate clustered heatmap and capture output
        ph <- pheatmap(expression_matrix_clean,
                    main = paste("DEG Patterns:", timepoint, "Treatments"),
                    color = heatmap_style$color_palette,
                    breaks = seq(-cap_value, cap_value, length.out = 51),
                    fontsize_row = 6,
                    fontsize_col = 9,
                    angle_col = heatmap_style$angle_col,
                    cluster_rows = TRUE,
                    cluster_cols = heatmap_style$cluster_cols,
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D2",
                    silent = TRUE)

        # Dynamic cluster number determination
        max_possible_clusters <- min(10, nrow(expression_matrix_clean)-1)
        clusters <- tryCatch({
            cutree(ph$tree_row, k = max_possible_clusters)
        }, error = function(e) {
            message("Cluster error: ", e$message)
            rep(1, nrow(expression_matrix_clean)) # Fallback to single cluster
        })

        # Functional annotation only for clusters with >5 genes
        if(length(unique(clusters)) > 1) {
            heatmap_fn <- file.path(heatmap_dir, paste0("DEG_heatmap_", timepoint, "_annotated.pdf"))

            cluster_terms <- lapply(unique(clusters), function(cl) {
                genes <- names(clusters[clusters == cl])
                if(length(genes) > 5) {
                gost_res <- gost(genes, organism = "athaliana", 
                                sources = c("GO:BP", "KEGG", "REAC"))
                if(!is.null(gost_res$result)) {
                    return(data.frame(
                    cluster = cl,
                    term = gost_res$result$term_name[1],
                    p = gost_res$result$p_value[1]
                    ))
                }
                }
                return(NULL)
            })

            cluster_df <- do.call(rbind, cluster_terms)
            if(!is.null(cluster_df)) {
                row_annot <- data.frame(
                    Function = sapply(clusters, function(x) {
                        cluster_df$term[cluster_df$cluster == x][1] %||% "Uncharacterized"
                    })
                )
                rownames(row_annot) <- names(clusters)
                
                # Generate annotated heatmap
                pdf(heatmap_fn, width = 14, height = 250)
                pheatmap(expression_matrix_clean,
                        annotation_row = row_annot,
                        main = paste("Functionally Annotated DEG Patterns:", timepoint),
                        show_rownames = heatmap_style$show_rownames,
                        gaps_row = cumsum(rle(clusters)$lengths)[-1],
                        cutree_rows = max_possible_clusters,
                        cluster_cols = heatmap_style$cluster_cols,
                        color = heatmap_style$color_palette,
                        border_color = heatmap_style$border_color,
                        cellheight = 10,
                        cellwidth = 20)
                dev.off()
                return()
            }
        }

        # Fallback: Generate basic heatmap if clustering fails
        pdf(heatmap_fn, width = 14, height = 24)
        pheatmap(expression_matrix_clean,
                main = paste("DEG Patterns:", timepoint, "Treatments"),
                show_rownames = heatmap_style$show_rownames,
                gaps_row = cumsum(rle(clusters)$lengths)[-1],
                cutree_rows = max_possible_clusters,
                color = heatmap_style$color_palette,
                cluster_cols = heatmap_style$cluster_cols,
                border_color = heatmap_style$border_color,
                cellheight = 10,
                cellwidth = 20)
        dev.off()
    }
}  

# Generate combined annotation report
if(exists("gene_annotations")) {
  write.csv(gene_annotations, file.path(base_output_dir, "combined_gene_annotations.csv"))
}
