library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(UpSetR)

# Set a directory for saving plots
plot_dir <- "./deseq2_output"  # Change this to your desired directory
dir.create(plot_dir, showWarnings = FALSE)  # Create the directory if it doesn't exist

# Function to print with timestamp and elapsed time
print_timed <- function(message, start_time = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (!is.null(start_time)) {
    elapsed <- Sys.time() - start_time
    message <- paste0(message, " (Completed in ", format(elapsed, digits = 3), ")")
  }
  cat(timestamp, ": ", message, "\n", sep = "")
}

# 1. Load data, skipping header lines and setting row names
counts <- read.table("./data/quant/clean_counts.txt", header = FALSE, skip = 2, row.names = 1, stringsAsFactors = FALSE, sep = "\t")

# 2. Read column names separately
column_names <- read.table("./data/quant/clean_counts.txt", header = FALSE, skip = 1, nrows = 1, stringsAsFactors = FALSE, sep = "\t")

# 3. Extract actual sample names, and remove bam from the name
sample_names <- gsub(".*/", "", column_names[1, -1])
sample_names <- gsub("\\.bam$", "", sample_names)

# 4. Assign sample names to counts
colnames(counts) <- sample_names

# 5. Load metadata
metadata <- read.csv("./metadata/metadata.csv", stringsAsFactors = FALSE)

# 6. Rename the metadata sample column
colnames(metadata)[colnames(metadata) == "sample_name"] <- "sample"

# 7. Update the metadata to drop the letter from HopN1a and HopB1a and the j from HopAB1j
metadata$treatment <- gsub("HopN1a", "HopN1", metadata$treatment)
metadata$treatment <- gsub("HopB1a", "HopB1", metadata$treatment)
metadata$treatment <- gsub("HopAB1j", "HopAB1", metadata$treatment)

# 8. Ensure 'treatment' and 'timepoint' are factors
metadata$treatment <- factor(metadata$treatment)
metadata$timepoint <- factor(metadata$timepoint)

# 9. Set reference level for hierarchical controls
metadata$treatment <- relevel(metadata$treatment, ref = "MgSO4")

# 10. Ensure metadata sample is character
metadata$sample <- as.character(metadata$sample)

# 12. Ensure that order is consistent and remove any rows that were dropped due to anti_join
metadata <- metadata[metadata$sample %in% colnames(counts),]

# 13. Make sure it's in the correct order
metadata <- metadata[match(colnames(counts), metadata$sample),]

# 14. Verify sample order and subset counts to match metadata
counts <- counts[, metadata$sample]

# 15. Print out the sample columns and the columns of count, before trying to subset
print("Metadata Sample Names:")
print(metadata$sample)
print("Count Matrix Column Names:")
print(colnames(counts))

# 16. Create DESeq2 object with interaction design
dds <- tryCatch({
    DESeqDataSetFromMatrix(countData = counts,
                           colData = metadata,
                           design = ~ treatment + timepoint + treatment:timepoint)
}, error = function(e) {
    stop("Error creating DESeq2 object: ", e$message)
})

# 17. Run DESeq analysis
start_time <- Sys.time()
print_timed("Running DESeq analysis...")
print_timed("  - Estimating size factors...")
dds <- estimateSizeFactors(dds)
print_timed("  - Estimating dispersion...")
dds <- estimateDispersions(dds)
print_timed("  - Fitting GLM and testing...")
dds <- nbinomWaldTest(dds)
print_timed("DESeq analysis complete", start_time)

# 18. Perform comparisons
print("Performing comparisons...")

    # List of active treatments (excluding D36E_EV and MgSO4)
    active_treatments <- levels(metadata$treatment)
    active_treatments <- active_treatments[! active_treatments %in% c("D36E_EV", "MgSO4")]

    # (a) Compare D36E::EV vs MgSO4 (global control comparison)
    print("  - Comparing D36E::EV vs MgSO4...")
    res_d36e_vs_mgso4 <- results(dds, contrast = c("treatment", "D36E_EV", "MgSO4"))
    write.csv(as.data.frame(res_d36e_vs_mgso4), file = file.path(plot_dir, "D36E_EV_vs_MgSO4.csv"))

    # Save significant genes for D36E::EV vs MgSO4
    sig_genes_d36e_vs_mgso4 <- subset(as.data.frame(res_d36e_vs_mgso4), padj < 0.05 & abs(log2FoldChange) > 1)
    write.csv(sig_genes_d36e_vs_mgso4, file = file.path(plot_dir, "Significant_Genes_D36E_EV_vs_MgSO4.csv"))

    # (b) Compare MgSO4 vs D36E::EV to filter out background noise
    print("  - Filtering out background noise from MgSO4...")
    res_mgso4 <- subset(as.data.frame(res_d36e_vs_mgso4), padj < 0.05 & abs(log2FoldChange) > 1)

    # Get the list of genes significantly expressed in MgSO4
    background_genes <- rownames(res_mgso4)

    # Filter out background genes from D36E::EV results
    filtered_res_d36e <- res_d36e_vs_mgso4[!rownames(res_d36e_vs_mgso4) %in% background_genes, ]
    write.csv(filtered_res_d36e, file = file.path(plot_dir, "Filtered_D36E_EV_vs_MgSO4.csv"))

    # Save filtered significant genes for D36E::EV
    filtered_sig_genes_d36e <- subset(filtered_res_d36e, padj < 0.05 & abs(log2FoldChange) > 1)
    write.csv(filtered_sig_genes_d36e, file = file.path(plot_dir, "Filtered_Significant_Genes_D36E_EV_vs_MgSO4.csv"))

    # (c) Compare each active treatment against filtered D36E::EV
    for (treatment in active_treatments) {
        print(paste("  - Comparing", treatment, "vs filtered D36E::EV..."))
        res <- results(dds, contrast = c("treatment", treatment, "D36E_EV"))
        
        # Remove background noise genes from active treatment results
        res_filtered <- res[!rownames(res) %in% background_genes, ]
        write.csv(as.data.frame(res_filtered), file = file.path(plot_dir, paste0(treatment, "_vs_Filtered_D36E_EV.csv")))
        
        # Save significant genes
        sig_genes <- subset(as.data.frame(res_filtered), padj < 0.05 & abs(log2FoldChange) > 1)
        write.csv(sig_genes, file = file.path(plot_dir, paste0("Significant_Genes_", treatment, "_vs_Filtered_D36E_EV.csv")))
        
        # Generate MA plot
        png(file.path(plot_dir, paste0("MA_plot_", treatment, "_vs_Filtered_D36E_EV.png")))
        plotMA(res_filtered, main = paste("MA Plot:", treatment, "vs Filtered D36E::EV"))
        dev.off()
        
        # Generate Volcano plot
        res_df <- as.data.frame(res_filtered)
        res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
        
        png(file.path(plot_dir, paste0("Volcano_plot_", treatment, "_vs_Filtered_D36E_EV.png")))
        print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
            geom_point() +
            theme_minimal() +
            labs(title = paste("Volcano Plot:", treatment, "vs Filtered D36E::EV"),
                x = "Log2 Fold Change",
                y = "-Log10 Adjusted P-value"))
        dev.off()
    }
print("Comparisons complete.")

# 19. Compare HopB1 vs HopAB1, excluding shared genes with D36E::EV
print("Comparing HopB1 vs HopAB1, excluding shared genes with D36E::EV...")

    # Get significant genes expressed in D36E::EV
    res_d36e_vs_mgso4 <- results(dds, contrast = c("treatment", "D36E_EV", "MgSO4"))
    sig_genes_d36e <- subset(as.data.frame(res_d36e_vs_mgso4), padj < 0.05 & abs(log2FoldChange) > 1)
    background_genes_d36e <- rownames(sig_genes_d36e)

    # Perform comparison between HopB1 and HopAB1
    res_hopb1_vs_hopab1 <- results(dds, contrast = c("treatment", "HopB1", "HopAB1"))

    # Exclude genes shared with D36E::EV
    filtered_res_hopb1_vs_hopab1 <- res_hopb1_vs_hopab1[!rownames(res_hopb1_vs_hopab1) %in% background_genes_d36e, ]
    write.csv(as.data.frame(filtered_res_hopb1_vs_hopab1), file = file.path(plot_dir, "Filtered_HopB1_vs_HopAB1.csv"))

    # Save significant genes for HopB1 vs HopAB1 after filtering
    filtered_sig_genes_hopb1_vs_hopab1 <- subset(filtered_res_hopb1_vs_hopab1, padj < 0.05 & abs(log2FoldChange) > 1)
    write.csv(filtered_sig_genes_hopb1_vs_hopab1, file = file.path(plot_dir, "Filtered_Significant_Genes_HopB1_vs_HopAB1.csv"))

    # Generate MA plot for filtered HopB1 vs HopAB1 comparison
    png(file.path(plot_dir, "MA_plot_Filtered_HopB1_vs_HopAB1.png"))
    plotMA(filtered_res_hopb1_vs_hopab1, main = "MA Plot: Filtered HopB1 vs HopAB1")
    dev.off()

    # Generate Volcano plot for filtered HopB1 vs HopAB1 comparison
    res_df <- as.data.frame(filtered_res_hopb1_vs_hopab1)
    res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")

    png(file.path(plot_dir, "Volcano_plot_Filtered_HopB1_vs_HopAB1.png"))
    print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
        geom_point() +
        theme_minimal() +
        labs(title = "Volcano Plot: Filtered HopB1 vs HopAB1",
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"))
    dev.off()

    print("Comparison between HopB1 and HopAB1 complete.")

# #20 Generate Enhanced Volcano Plots for filtered results
# enhanced_volcano_dir = "./enhanced_volcano_plots"
# for (treatment in active_treatments) {
#     res_filtered <- read.csv(file.path(plot_dir, paste0(treatment, "_vs_Filtered_D36E_EV.csv")), row.names = 1)
    
#     png(file.path(plot_dir, paste0("Enhanced_Volcano_", treatment, ".png")))
#     tryCatch({
#         EnhancedVolcano(res_filtered,
#             lab = rownames(res_filtered),
#             x = 'log2FoldChange',
#             y = 'padj',
#             title = 'Enhanced Volcano Plot',
#             pCutoff = 0.05,
#             FCcutoff = 1,
#             pointSize = 3.0,
#             labSize = 4.0,
#             colAlpha = 0.7)
#     }, error = function(e) {
#         warning("Error generating Enhanced Volcano plot: ", e$message)
#     })
#     dev.off()
# }

# #21 EFFECTOR-SPECIFIC FINGERPRINT ANALYSIS
# # To identify shared and unique DEGs between treatments.
#     # Load significant genes for each treatment
#     sig_genes_list <- lapply(active_treatments, function(treatment) {
#         read.csv(file.path(plot_dir, paste0("Significant_Genes_", treatment, "_vs_Filtered_D36E_EV.csv")), row.names = 1)$gene
#     })

#     names(sig_genes_list) <- active_treatments

#     # Generate UpSet plot to visualize overlaps
#     pdf(file.path(plot_dir, "Effector_Fingerprint_UpSet.pdf"))
#     upset(fromList(sig_genes_list), 
#         main.bar.color = "blue", 
#         sets.bar.color = "red",
#         order.by = "freq",
#         text.scale = c(1.5, 1.5, 1.2))
#     dev.off()


# #22 TIME-RESOLVED ANALYSIS
# # Split DEG analysis by timepoint (1h vs. 8h)
# # Perform timepoint-specific comparisons
# for (timepoint in levels(metadata$timepoint)) {
#     for (treatment in active_treatments) {
#         res_timepoint <- results(dds, name = paste0("treatment", treatment, ".timepoint", timepoint))
        
#         # Save results and significant genes
#         write.csv(as.data.frame(res_timepoint), file = file.path(plot_dir, paste0(treatment, "_", timepoint, "_vs_Filtered_D36E_EV.csv")))
        
#         sig_genes_timepoint <- subset(as.data.frame(res_timepoint), padj < 0.05 & abs(log2FoldChange) > 1)
#         write.csv(sig_genes_timepoint, file = file.path(plot_dir, paste0("Significant_Genes_", treatment, "_", timepoint, "_vs_Filtered_D36E_EV.csv")))
#     }
# }

# #23 CO-EXPRESSION NETWORK ANALYSIS
# # To identify effector-specific gene modules.
# # Integrate WGCNA to identify gene modules correlated with specific effectors:
#     # Prepare data for WGCNA
#     expr_data <- counts[rownames(counts) %in% unique(unlist(sig_genes_list)), ]
#     trait_data <- metadata[, c("sample", "treatment", "timepoint")]

#     # Run WGCNA pipeline
#     net <- blockwiseModules(expr_data,
#                             power = 6,
#                             TOMType = "unsigned",
#                             minModuleSize = 30,
#                             reassignThreshold = 0,
#                             mergeCutHeight = 0.25,
#                             numericLabels = TRUE,
#                             pamRespectsDendro = FALSE,
#                             saveTOMs = TRUE,
#                             saveTOMFileBase = "WGCNA_TOM",
#                             verbose = 3)

#     # Relate modules to traits
#     moduleTraitCor <- cor(net$MEs, trait_data$treatment)
#     write.csv(moduleTraitCor, file.path(plot_dir, "Module_Trait_Correlations.csv"))

# #24 VALIDATE KEY FINGERPRINT GENES
# # Prioritize highly significant genes and cross-reference with known effector targets
#     # Filter top candidates for validation
#     top_candidates <- subset(res_filtered, padj < 0.01 & abs(log2FoldChange) > 2)
#     write.csv(top_candidates, file.path(plot_dir, "Top_Candidate_Genes.csv"))

#     # Cross-reference with known targets (e.g., JAZ proteins or PR genes)
#     known_targets <- c("JAZ1", "PR1", "PDF1.2")
#     validated_candidates <- top_candidates[top_candidates$gene %in% known_targets, ]
#     write.csv(validated_candidates, file.path(plot_dir, "Validated_Candidates.csv"))
