library(pheatmap)
library(gprofiler2)
library(reshape2)
library(gridExtra)
library(purrr)

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

    ##################################################
    ############ DEG HEATMAP GENERATION ##############
    ##################################################

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







    ##################################################
    ############## ENHANCED GO ANALYSIS ##############
    ##################################################

        # Create GO output directory
        go_dir <- file.path(base_output_dir, "go_enrichment", timepoint)
        dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Process each treatment with robust error handling
        for(file in current_files) {
            treatment <- gsub("SigGenes_|\\.csv", "", basename(file))
            treatment_clean <- gsub("_vs_D36E_EV", "", treatment)
            
            # Read and validate DEG data
            deg_data <- tryCatch({
                raw <- read.csv(file)
                colnames(raw)[1] <- "gene_id"
                
                # Validate critical columns
                required_cols <- c("gene_id", "log2FoldChange", "padj")
                if(!all(required_cols %in% colnames(raw))) {
                    stop("Missing required columns in ", file)
                }
                
                # Clean and filter genes
                raw$gene_id <- gsub("\\.\\d+$", "", raw$gene_id)
                valid <- raw[raw$padj < CRITICAL_VALUE & 
                            abs(raw$log2FoldChange) > lfc_cutoff, ]
                
                if(nrow(valid) == 0) {
                    message("No DEGs for ", treatment_clean)
                    next
                }
                
                valid$Direction <- ifelse(valid$log2FoldChange > 0, "Up", "Down")
                valid
            }, error = function(e) {
                message("Error: ", e$message)
                return(NULL)
            })
            
            if(is.null(deg_data)) next
            
            # GO analysis with API retries
            perform_go <- function(genes) {
                attempt <- 1
                while(attempt <= 3) {
                    result <- tryCatch({
                        gprofiler2::gost(
                            genes,
                            organism = "athaliana",
                            sources = c("GO:BP", "KEGG", "REAC"),
                            significant = TRUE,
                            user_threshold = 0.05,
                            correction_method = "g_SCS",
                            evcodes = TRUE  # Critical for intersection data [3][5]
                        )$result
                    }, error = function(e) {
                        message("API attempt ", attempt, " failed: ", e$message)
                        Sys.sleep(5)
                        return(NULL)
                    })
                    
                    if(!is.null(result)) break
                    attempt <- attempt + 1
                }
                return(result)
            }
            
            # Process both directions with column validation
            results <- list()
            for(direction in c("Up", "Down")) {
                genes <- deg_data$gene_id[deg_data$Direction == direction]
                
                # Adaptive thresholds based on treatment type
                min_genes <- if(grepl("HopN1", treatment_clean)) {
                    3
                } else if(timepoint == "1h") {
                    5
                } else {
                    10
                }
                
                if(length(genes) < min_genes) {
                    message(treatment_clean, " ", direction, ": ", length(genes), " genes")
                    next
                }
                
                # Get GO results with validation
                go_res <- perform_go(genes)
                if(is.null(go_res) || nrow(go_res) == 0) next
                
                # Check for required columns
                required_cols <- c("term_id", "term_name", "p_value", "source", 
                                "intersection_size", "query_size", "intersection")
                missing_cols <- setdiff(required_cols, colnames(go_res))
                if(length(missing_cols) > 0) {
                    message("Missing columns in GO results: ", paste(missing_cols, collapse=", "))
                    next
                }
                
                # Filter and format results
                clean_terms <- go_res[!grepl("^GO:0008150|^GO:0005575|^GO:0003674", go_res$term_id), ]
                if(nrow(clean_terms) == 0) next
                
                # Select top terms per source
                sources <- unique(clean_terms$source)
                top_terms <- do.call(rbind, lapply(sources, function(s) {
                    subset <- clean_terms[clean_terms$source == s, ]
                    subset[order(subset$p_value), ][1:min(15, nrow(subset)), ]
                }))
                
                # Add metadata with validation
                top_terms <- cbind(
                    top_terms,
                    data.frame(
                        Treatment = treatment_clean,
                        Timepoint = timepoint,
                        Direction = direction,
                        GeneRatio = paste0(top_terms$intersection_size, "/", top_terms$query_size),
                        stringsAsFactors = FALSE
                    )
                )
                
                results[[direction]] <- top_terms
            }
            
            # Safe CSV writing with column validation
            output_file <- file.path(go_dir, paste0(treatment_clean, "_GO.csv"))
            if(length(results) > 0) {
                final <- do.call(rbind, results)
                # Validate output columns exist
                output_cols <- c("term_id", "term_name", "p_value", "source",
                                "GeneRatio", "intersection", "Treatment",
                                "Timepoint", "Direction")
                existing_cols <- intersect(output_cols, colnames(final))
                
                if(length(existing_cols) > 0) {
                    write.csv(final[, existing_cols], output_file, row.names = FALSE)
                } else {
                    write.csv(data.frame(Note = "No valid columns found"), 
                            output_file, row.names = FALSE)
                }
            } else {
                write.csv(data.frame(Note = "No significant terms found"), 
                        output_file, row.names = FALSE)
            }
        }











    ##################################################
    ############# CREATE GO TERM HEATMAP #############
    ##################################################

    # Check if GO results exist for this timepoint
    go_dir <- file.path(base_output_dir, "go_enrichment", timepoint)
    go_files <- list.files(go_dir, pattern = "_GO.csv", full.names = TRUE)
    
    if(length(go_files) > 0) {
        # Read and process GO results with robust error handling
        go_data <- do.call(rbind, lapply(go_files, function(f) {
            df <- tryCatch({
                temp <- read.csv(f)
                # Create unique biological identifiers with validation
                temp$unique_id <- ifelse(
                    is.na(temp$term_id) | is.na(temp$term_name),
                    NA,
                    paste0(temp$term_id, ": ", temp$term_name)
                )
                # Filter invalid entries and root terms
                temp <- temp[!is.na(temp$unique_id) & 
                            !grepl("^GO:0008150|^GO:0005575|^GO:0003674", temp$term_id), ]
                
                if(nrow(temp) == 0) return(NULL)
                
                # Calculate log p-values with biological constraints
                temp$log_p <- -log10(temp$p_value + 1e-100)
                # temp[is.infinite(temp$log_p), "log_p"] <- 20  # Cap extreme values
                temp[, c("unique_id", "Treatment", "log_p")]
            }, error = function(e) {
                message("Error reading ", f, ": ", e$message)
                return(NULL)
            })
        }))
        
        # Create heatmap matrix if valid biological data exists
        if(!is.null(go_data) && nrow(go_data) > 0) {
            # Create matrix with enforced unique identifiers
            heatmap_matrix <- reshape2::dcast(go_data, unique_id ~ Treatment, 
                                            value.var = "log_p", fun.aggregate = max)
            
            # Remove invalid rows using biological thresholds
            valid_rows <- apply(heatmap_matrix[,-1], 1, function(x) {
                sum(!is.na(x)) >= 1 &&         # Minimum biological observations
                var(x, na.rm = TRUE) > 0.1 &&  # Meaningful biological variation
                all(is.finite(x))              # Remove non-finite values
            })
            heatmap_matrix <- heatmap_matrix[valid_rows & !is.na(heatmap_matrix$unique_id), ]
            
            if(nrow(heatmap_matrix) > 2 && ncol(heatmap_matrix) > 1) {
                # Set validated row names with biological identifiers
                rownames(heatmap_matrix) <- make.unique(
                    as.character(heatmap_matrix$unique_id), 
                    sep = "_"
                )
                heatmap_matrix <- as.matrix(heatmap_matrix[,-1])
                
                # Apply biological value constraints
                heatmap_matrix[is.na(heatmap_matrix)] <- 0
                heatmap_matrix[heatmap_matrix > 20] <- 20  # Cap at p < 1e-20
                
                # Generate heatmap with failsafe parameters
                go_heatmap_fn <- file.path(go_dir, paste0("GO_term_heatmap_", timepoint, ".pdf"))
                pdf(go_heatmap_fn, width = 10, height = max(6, nrow(heatmap_matrix) * 0.3))
                
                pheatmap(heatmap_matrix,
                        main = paste("GO Term Enrichment:", timepoint),
                        breaks = seq(0, 20, length.out = 51),
                        fontsize_row = 7,
                        fontsize_col = 9,
                        angle_col = 45,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        cellwidth = 30,
                        cellheight = 15,
                        clustering_method = "ward.D2",
                        border_color = NA,
                        show_rownames = TRUE,
                        color = colorRampPalette(c("white", "blue"))(50))
                dev.off()
            } else {
                message("Insufficient valid biological terms for clustering in ", timepoint)
            }
        }
    }












    
}

















##################################################
############ REGENERATE DEG HEATMAPS #############
##################################################


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
    expression_matrix_clean <- expression_matrix_clean[rowSums(!is.na(expression_matrix)) >= 1, ]

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