library(gprofiler2)
library(ggplot2)
library(pheatmap)

# Define analysis parameters
CRITICAL_VALUE <- 0.05
padj_cutoff <- CRITICAL_VALUE
lfc_cutoff <- 1 

# Specify the base output directory
base_output_dir <- "./results"
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# List of CSV files
deg_files <- c("./deseq2_output/joint_timepoints/Significant_Genes_HopAB1_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/joint_timepoints/Significant_Genes_HopB1_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/joint_timepoints/Significant_Genes_HopN1_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopAB1_1h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopAB1_8h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopB1_1h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopB1_8h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopN1_1h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/time_resolved_analysis/Significant_Genes_HopN1_8h_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/within_treatment_time_comparisons/Significant_Genes_HopAB1_8h_vs_1h.csv",
               "./deseq2_output/within_treatment_time_comparisons/Significant_Genes_HopB1_8h_vs_1h.csv",
               "./deseq2_output/within_treatment_time_comparisons/Significant_Genes_HopN1_8h_vs_1h.csv")

for (file in deg_files) {
    # print(paste("=== Processing file:", file, "==="))
    # Create output subdirectory mirroring input structure
    output_subdir <- file.path(base_output_dir, dirname(gsub("^\\./deseq2_output/", "", file)))
    dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
    
    # Read and filter DEG data
    deg_data <- read.csv(file, header = TRUE)
    colnames(deg_data)[1] <- "gene_id"
    # print(paste("Columns detected:", paste(colnames(deg_data), collapse = ", ")))
    
    sig_genes <- deg_data[deg_data$padj < padj_cutoff & 
                        abs(deg_data$log2FoldChange) > lfc_cutoff, ]
    print(paste(file, paste("|", nrow(sig_genes), "significant genes")))

    # Perform functional enrichment
    gost_results <- gost(query = sig_genes$gene_id,
                        organism = "athaliana",
                        sources = c("GO", "KEGG", "REAC"),
                        evcodes = TRUE,
                        correction_method = "fdr"
    )

    # Save full results
    output_prefix <- file.path(output_subdir, tools::file_path_sans_ext(basename(file)))

    # Convert and save gprofiler results
    gost_df <- as.data.frame(gost_results$result)
    list_cols <- sapply(gost_df, is.list)
    gost_df[list_cols] <- lapply(gost_df[list_cols], function(x) sapply(x, paste, collapse = ","))
    write.csv(gost_df, paste0(output_prefix, "_gprofiler_enrichment.csv"), row.names = FALSE)

    # Generate Manhattan plot
    pdf(paste0(output_prefix, "_enrichment_plot.pdf"))
    plot_title <- paste("Enrichment Analysis:", tools::file_path_sans_ext(basename(file)))
    print(gostplot(gost_results, capped = FALSE, interactive = FALSE) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    dev.off()

    # Generate Heatmap
    if(nrow(sig_genes) > 1 && !is.null(gost_results$result)) {
        tryCatch({
            # Extract terms and genes
            terms <- gost_results$result$term_id
            genes <- unique(unlist(strsplit(gost_results$result$intersection, ",")))
            
            # Create matrix
            heatmap_matrix <- matrix(0, nrow = length(genes), ncol = length(terms))
            rownames(heatmap_matrix) <- genes
            colnames(heatmap_matrix) <- terms
            
            # Fill matrix
            for(i in 1:nrow(gost_results$result)) {
                term_genes <- unlist(strsplit(gost_results$result$intersection[i], ","))
                heatmap_matrix[term_genes, gost_results$result$term_id[i]] <- 1
            }
            
            # Generate heatmap
            heatmap_fn <- paste0(output_prefix, "_enrichment_heatmap.pdf")
            pdf(heatmap_fn, width = 12, height = 8)  # Adjust size as needed
            
            pheatmap(heatmap_matrix,
                    main = paste("Enrichment Heatmap:", tools::file_path_sans_ext(basename(file))),
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    fontsize_row = 6,
                    fontsize_col = 8,
                    angle_col = 45,
                    treeheight_row = 20,
                    treeheight_col = 20,
                    color = colorRampPalette(c("white", "firebrick3"))(2))
            
            dev.off()
            # print(paste("Successfully generated enrichment heatmap:", heatmap_fn))
        }, error = function(e) {
            message(paste("Error generating enrichment heatmap for", file))
            message(e)
        })
    } else {
        print(paste("Skipping heatmap - insufficient genes or enrichment results"))
    }


    # Additional annotation
    gene_annotations <- gconvert(
        sig_genes$gene_id,
        organism = "athaliana",
        target = "ENSG",
        numeric_ns = "ENTREZGENE_ACC"
    )
    write.csv(gene_annotations, paste0(output_prefix, "_gene_annotations.csv"))

    # Save cache
    saveRDS(gost_results, paste0(output_prefix, "_enrichment_cache.rds"))
}








COLOR_PALETTE = c("white", "dark blue")
COLOR_RESOLUTION = 5

# Create heatmap directory
heatmap_dir <- file.path(base_output_dir, "heatmaps")
dir.create(heatmap_dir, showWarnings = FALSE)



extract_enriched_terms <- function(file) {
    csv_file <- file.path("./results", 
                        dirname(gsub("^\\./deseq2_output/", "", file)), 
                        paste0(tools::file_path_sans_ext(basename(file)), "_gprofiler_enrichment.csv"))

    if (file.exists(csv_file)) {
        gost_results <- read.csv(csv_file)
        if (nrow(gost_results) > 0) {
            treatment <- sub(".*Significant_Genes_(.+)_vs.*", "\\1", basename(file))
            return(data.frame(term = gost_results$term_name, 
                            p_value = as.numeric(gost_results$p_value),
                            treatment = treatment))
        }
    }
    return(NULL)
}

filtered_deg_files <- deg_files[!grepl("joint_timepoints", deg_files)]

all_enriched_terms <- do.call(rbind, lapply(filtered_deg_files, extract_enriched_terms))

if (!is.null(all_enriched_terms) && nrow(all_enriched_terms) > 0) {
    # Ensure p_value is numeric and remove any NA values
    all_enriched_terms$p_value <- as.numeric(all_enriched_terms$p_value)
    all_enriched_terms <- all_enriched_terms[!is.na(all_enriched_terms$p_value), ]
    
    # Create a matrix of -log10(p-values), using min as the aggregation function
    heatmap_data <- dcast(all_enriched_terms, term ~ treatment, value.var = "p_value", 
                          fun.aggregate = min, fill = 1)
    rownames(heatmap_data) <- heatmap_data$term
    heatmap_data$term <- NULL
    heatmap_data <- -log10(heatmap_data)
    
    # Sort columns alphabetically
    heatmap_data <- heatmap_data[, sort(colnames(heatmap_data))]
    
    # Determine optimal number of clusters using Elbow Method
    library(cluster)

    dist_matrix <- dist(heatmap_data)  # Compute distance matrix
    hclust_result <- hclust(dist_matrix)  # Perform hierarchical clustering

    # Calculate within-cluster sum of squares for different numbers of clusters
    wss <- sapply(1:80, function(k) {
        clusters <- cutree(hclust_result, k = k)
        sum(sapply(unique(clusters), function(cluster) {
            cluster_points <- heatmap_data[clusters == cluster, ]
            sum((cluster_points - colMeans(cluster_points))^2)
        }))
    })

    # Plot WSS to find elbow point
    elbow_plot_fn <- file.path(heatmap_dir, "elbow_plot.pdf")
    pdf(elbow_plot_fn)
    plot(1:80, wss, type = "b", pch = 19, frame = FALSE,
         xlab = "Number of Clusters",
         ylab = "Within-Cluster Sum of Squares",
         main = "Elbow Method for Determining Optimal Clusters")
    dev.off()
    print(paste("Generated elbow plot:", elbow_plot_fn))

    # Set number of clusters based on visual inspection or statistical determination
    num_clusters <- 40  # Adjust based on elbow plot or other methods

    # Generate comprehensive heatmap and save clustering information
    comprehensive_heatmap_fn <- file.path(heatmap_dir, "comprehensive_enrichment_heatmap_time_resolved.png")

    # Set PNG output
    png(comprehensive_heatmap_fn, width = length(filtered_deg_files) * 450, height = max(3000, nrow(heatmap_data) * 15), res = 300)

    pheatmap_result <- pheatmap(
        heatmap_data,
        main = "Enrichment Heatmap by GO Term Across Time-Resolved Treatments (-log10(p-value))",
        show_rownames = TRUE,
        show_colnames = TRUE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        fontsize_row = 6,
        fontsize_col = 8,
        angle_col = 45,
        treeheight_row = 20,
        treeheight_col = 0,
        cellheight = 3.5,  # Adjust row height for better readability
        cellwidth = 25,   # Adjust column width for better spacing
        color = colorRampPalette(COLOR_PALETTE)(COLOR_RESOLUTION),
        legend_labels = "-log10(p-value)",
        legend = TRUE
    )

    dev.off()
    print(paste("Generated comprehensive enrichment heatmap:", comprehensive_heatmap_fn))

        
        # Load the specified terms and their clusters
    specified_terms <- c(
                            "response to stimulus",
                            "response to chemical",
                            "response to stress",
                            "bounding membrane of organelle",
                            "chloroplast envelope",
                            "purine nucleotide Binding",
                            "burine ribonucleotide binding",
                            "Carbohydrate derivative binding",
                            "ribonucleotide binding",
                            "adenyl ribonucleotide binding",
                            "kinase activity",
                            "intracellular anatomical structure",
                            "heterocyclic compound binding",
                            "Nucleoside phosphate binding",
                            "Anion binding",
                            "protein kinase activity",
                            "protein serine kinase activity",
                            "phosphotransferase activity",
                            "alcohol group as acceptor",
                            "protein serine/threonine kinase activity",
                            "ADP binding",
                            "immune response",
                            "immune system process",
                            "defense response to symbiont",
                            "innate immune response",
                            "protein modification process",
                            "metabolic process",
                            "protein phosphorylation",
                            "photosynthesis, light reaction",
                            "photosystem I",
                            "photosystem II",
                            "Plastid",
                            "Cytoplasm",
                            "organelle quter membrane",
                            "outer membrane",

                            "Response to endogenous stimulus",
                            "Response to hormone",
                            "Response to lipid",
                            "Response to wounding",
                            "defense response to other organism",
                            "signal transduction",
                            "cell communication",
                            "response to bacterium",
                            "signaling",
                            "cellular response to stress",
                            "Defense response to bacterium",
                            "regulation of cellular process",
                            "biological regulation",
                            "regulation of biological process",
                            "cellular process",
                            "photosynthesis",
                            "photosystem",
                            "plastid membrane",
                            "thylakoid",
                            "chloroplast thylakoid",
                            "Plastid envelope",
                            "photosynthetic membrane",
                            "thylakoid membrane",
                            "chloroplast thylakoid membrane",
                            "plastid thylakoid membrane",
                            "Cellular anatomical structure",
                            "cell periphery",
                            "plasma membrane",
                            "catalytic activity",
                            "carbohydrate biosynthetic process",
                            "membrane",
                            "Response to abiotic stimulus",
                            "response to oxygen-containing compound",
                            "Cellular response to chemical stimulus",
                            "Response to external stimulus",
                            "response to biotic stimulus",
                            "Biological process involved in interspecies interaction between organisms",
                            "response to other organism",
                            "cellular response to Hypoxia",
                            "cellular response to decreased oxygen levels",
                            "cellular response to oxygen levels",
                            "response to hypoxia",
                            "response to decreased oxygen levels",
                            "cellular response to stimulus",
                            "defense response"
    )

    # Filter heatmap data to include only specified terms
    filtered_heatmap_data <- heatmap_data[rownames(heatmap_data) %in% specified_terms, ]

    # Generate cluster-specific heatmap for the filtered terms
    cluster_specific_heatmap_fn <- file.path(heatmap_dir, "cluster_specific_enrichment_heatmap.png")
    png(cluster_specific_heatmap_fn, width = length(filtered_deg_files) * 450, height = max(3000, nrow(filtered_heatmap_data) * 75), res = 300)

    # Define legend breaks and labels
    # legend_breaks <- seq(0, max(heatmap_data, na.rm = TRUE), length.out = 5)  # Dynamically set breaks
    # legend_labels <- c("-log10(p-value) = 0", "Medium", "High", "Very High", "Extremely High")

    pheatmap(
        filtered_heatmap_data,
        main = "Enrichment Heatmap for Selected GO Terms (-log10(p-value))",
        show_rownames = TRUE,
        show_colnames = TRUE,
        cluster_rows = TRUE,  # Disable clustering for rows (already filtered)
        cluster_cols = FALSE,
        fontsize_row = 8,
        fontsize_col = 8,
        angle_col = 45,
        cellheight = 12.5,  # Adjust row height for readability
        cellwidth = 25,   # Adjust column width for better spacing
        color = colorRampPalette(COLOR_PALETTE)(COLOR_RESOLUTION),
        # legend_breaks = legend_breaks
        # legend_labels = legend_labels
    )

    dev.off()
    print(paste("Generated cluster-specific enrichment heatmap:", cluster_specific_heatmap_fn))
    } else {
        print(paste("No terms found in Cluster", cluster_of_interest))
}


















# extract_enriched_terms <- function(file) {
#     csv_file <- file.path("./results", 
#                         dirname(gsub("^\\./deseq2_output/", "", file)), 
#                         paste0(tools::file_path_sans_ext(basename(file)), "_gprofiler_enrichment.csv"))

#     if (file.exists(csv_file)) {
#     gost_results <- read.csv(csv_file)
#     if (nrow(gost_results) > 0) {
#         treatment <- sub(".*Significant_Genes_(.+)_vs.*", "\\1", basename(file))
#         return(data.frame(term = gost_results$term_name, 
#                         p_value = as.numeric(gost_results$p_value),
#                         treatment = treatment))
#     }
#     }
#     return(NULL)
# }

# filtered_deg_files <- deg_files[!grepl("joint_timepoints", deg_files)]

# all_enriched_terms <- do.call(rbind, lapply(filtered_deg_files, extract_enriched_terms))

# if (!is.null(all_enriched_terms) && nrow(all_enriched_terms) > 0) {
#     # Ensure p_value is numeric and remove any NA values
#     all_enriched_terms$p_value <- as.numeric(all_enriched_terms$p_value)
#     all_enriched_terms <- all_enriched_terms[!is.na(all_enriched_terms$p_value), ]

#     # Create a matrix of -log10(p-values), using min as the aggregation function
#     heatmap_data <- dcast(all_enriched_terms, term ~ treatment, value.var = "p_value", 
#                         fun.aggregate = min, fill = 1)
#     rownames(heatmap_data) <- heatmap_data$term
#     heatmap_data$term <- NULL
#     heatmap_data <- -log10(heatmap_data)

#     # Sort columns alphabetically
#     heatmap_data <- heatmap_data[, sort(colnames(heatmap_data))]

#     # Generate comprehensive heatmap
#     heatmap_fn <- paste0(base_output_dir, "/comprehensive_enrichment_heatmap_time_resolved.pdf")
#     pdf(heatmap_fn, width = 20, height = 30)  # Adjust size as needed

#     # Generate heatmap and save clustering information
#     pheatmap_result <- pheatmap(
#     heatmap_data,
#     main = "Comprehensive Enrichment Heatmap Across Time-Resolved Treatments",
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     cluster_rows = TRUE,
#     cluster_cols = FALSE,
#     fontsize_row = 6,
#     fontsize_col = 8,
#     angle_col = 45,
#     treeheight_row = 20,
#     treeheight_col = 0,
#     color = colorRampPalette(c("white", "light grey", "dark grey", "black"))(2),
#     legend_labels = "-log10(p-value)",
#     legend = TRUE
#     )

#     # Extract row clusters
#     num_clusters <- 4  # Specify the number of clusters you want
#     row_clusters <- cutree(pheatmap_result$tree_row, k = num_clusters)

#     # Save cluster assignments to a data frame
#     cluster_assignments <- data.frame(
#     term = rownames(heatmap_data),
#     cluster = row_clusters
#     )

#     # Write cluster assignments to a CSV file for reference
#     write.csv(cluster_assignments, file.path(base_output_dir, "term_cluster_assignments.csv"), row.names = FALSE)

#     dev.off()
#     print(paste("Generated comprehensive enrichment heatmap for time-resolved analysis:", heatmap_fn))
# } else {
#   print("No enriched terms found across time-resolved treatments")
# }

























































extract_enriched_terms <- function(file) {
  csv_file <- file.path("./results", 
                        dirname(gsub("^\\./deseq2_output/", "", file)), 
                        paste0(tools::file_path_sans_ext(basename(file)), "_gprofiler_enrichment.csv"))
  
  if (file.exists(csv_file)) {
    gost_results <- read.csv(csv_file)
    if (nrow(gost_results) > 0) {
      treatment <- sub(".*Significant_Genes_(.+)_vs.*", "\\1", basename(file))
      return(data.frame(source = gost_results$source, 
                        p_value = as.numeric(gost_results$p_value),
                        treatment = treatment))
    }
  }
  return(NULL)
}

filtered_deg_files <- deg_files[!grepl("joint_timepoints", deg_files)]

all_enriched_terms <- do.call(rbind, lapply(filtered_deg_files, extract_enriched_terms))

if (!is.null(all_enriched_terms) && nrow(all_enriched_terms) > 0) {
  all_enriched_terms$p_value <- as.numeric(all_enriched_terms$p_value)
  all_enriched_terms <- all_enriched_terms[!is.na(all_enriched_terms$p_value), ]
  
  # Aggregate data by source
  aggregated_data <- aggregate(p_value ~ source + treatment, data = all_enriched_terms, FUN = min)
  
  # Create a matrix of -log10(p-values)
  heatmap_data <- dcast(aggregated_data, source ~ treatment, value.var = "p_value", fill = NA)
  rownames(heatmap_data) <- heatmap_data$source
  heatmap_data$source <- NULL
  heatmap_data[is.na(heatmap_data)] <- 1
  heatmap_data <- -log10(heatmap_data)
  
  # Sort columns alphabetically
  heatmap_data <- heatmap_data[, sort(colnames(heatmap_data))]
  
  # Generate comprehensive heatmap
  heatmap_fn <- file.path(heatmap_dir, "comprehensive_enrichment_heatmap_by_source.pdf")
  pdf(heatmap_fn, width = 12, height = 8)  # Adjusted size for fewer rows
  
  pheatmap(heatmap_data,
           main = "Enrichment Heatmap by GO Source Across Time-Resolved Treatments (-log10(p-value))",
           show_rownames = TRUE,
           show_colnames = TRUE,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           cellwidth = 30,
           cellheight = 30,
           color = colorRampPalette(COLOR_PALETTE)(COLOR_RESOLUTION),
           legend_labels = "-log10(p-value)",
           legend = TRUE)
  
  dev.off()
  print(paste("Generated comprehensive enrichment heatmap by GO source:", heatmap_fn))
}
