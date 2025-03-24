# Load required libraries (no Bioconductor dependencies)
library(gprofiler2)
library(ggplot2)

# Define analysis parameters
CRITICAL_VALUE <- 0.05
padj_cutoff <- CRITICAL_VALUE
lfc_cutoff <- 1 

# Specify the output directory
output_dir <- "./results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List of CSV files
deg_files <- c("./deseq2_output/Significant_Genes_HopAB1_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/Significant_Genes_HopB1_vs_Filtered_D36E_EV.csv",
               "./deseq2_output/Significant_Genes_HopN1_vs_Filtered_D36E_EV.csv")

for (file in deg_files) {
    print(file)
    # Read and filter DEG data
    deg_data <- read.csv(file, header = TRUE)
    colnames(deg_data)[1] <- "gene_id"
    print(colnames(deg_data))
    sig_genes <- deg_data[deg_data$padj < padj_cutoff & 
                        abs(deg_data$log2FoldChange) > lfc_cutoff, ]

    # Perform functional enrichment using g:Profiler web service
    gost_results <- gost(query = sig_genes$gene_id,
                        organism = "athaliana",  # Arabidopsis thaliana
                        sources = c("GO", "KEGG", "REAC"),
                        evcodes = TRUE,  # Include evidence codes
                        correction_method = "fdr"
    )

    # Save full results
    output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(file)))

    # Convert gost_results$result to a data frame and flatten list columns
    gost_df <- as.data.frame(gost_results$result)
    list_cols <- sapply(gost_df, is.list)
    gost_df[list_cols] <- lapply(gost_df[list_cols], function(x) sapply(x, paste, collapse = ","))

    # Write the flattened data frame to CSV
    write.csv(gost_df, paste0(output_prefix, "_gprofiler_enrichment.csv"), row.names = FALSE)

    # Generate Manhattan plot
    pdf(paste0(output_prefix, "_enrichment_plot.pdf"))
    print(gostplot(gost_results, capped = FALSE, interactive = FALSE))
    dev.off()

    # Additional annotation using built-in functions
    gene_annotations <- gconvert(
    sig_genes$gene_id,
    organism = "athaliana",
    target = "ENSG",
    numeric_ns = "ENTREZGENE_ACC"
    )
    write.csv(gene_annotations, paste0(output_prefix, "_gene_annotations.csv"))

    # Save annotation cache for later offline use
    saveRDS(gost_results, paste0(output_prefix, "_enrichment_cache.rds"))
}
