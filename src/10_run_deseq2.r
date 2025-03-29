library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(UpSetR)

# Constants
CRITICAL_VALUE <- 0.05
MIN_LOG2FOLD_CHANGE <- 1

# Set a directory for saving plots
working_dir <- "./deseq2_output"  # Change this to your desired directory
dir.create(working_dir, showWarnings = FALSE)  # Create the directory if it doesn't exist

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

# 7. Update the metadata to drop the letter a from HopN1a and HopB1a and the j from HopAB1j
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

    joint_output_dir <- file.path(working_dir, "/joint_timepoints")
    dir.create(joint_output_dir, showWarnings = FALSE)

    # List of active treatments (excluding D36E_EV and MgSO4)
    active_treatments <- levels(metadata$treatment)
    active_treatments <- active_treatments[! active_treatments %in% c("D36E_EV", "MgSO4")]

    # (a) Compare D36E::EV vs MgSO4 (global control comparison)
    print("  - Comparing D36E::EV vs MgSO4...")
    res_d36e_vs_mgso4 <- results(dds, contrast = c("treatment", "D36E_EV", "MgSO4"))
    write.csv(as.data.frame(res_d36e_vs_mgso4), file = file.path(joint_output_dir, "D36E_EV_vs_MgSO4.csv"))

    # Save significant genes for D36E::EV vs MgSO4
    sig_genes_d36e_vs_mgso4 <- subset(as.data.frame(res_d36e_vs_mgso4), padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
    write.csv(sig_genes_d36e_vs_mgso4, file = file.path(joint_output_dir, "Significant_Genes_D36E_EV_vs_MgSO4.csv"))

    # (b) Compare MgSO4 vs D36E::EV to filter out background noise
    print("  - Filtering out background noise from MgSO4...")
    res_mgso4 <- subset(as.data.frame(res_d36e_vs_mgso4), padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)

    # Get the list of genes significantly expressed in MgSO4
    background_genes <- rownames(res_mgso4)

    # Filter out background genes from D36E::EV results
    filtered_res_d36e <- res_d36e_vs_mgso4[!rownames(res_d36e_vs_mgso4) %in% background_genes, ]
    write.csv(filtered_res_d36e, file = file.path(joint_output_dir, "Filtered_D36E_EV_vs_MgSO4.csv"))

    # Save filtered significant genes for D36E::EV
    filtered_sig_genes_d36e <- subset(filtered_res_d36e, padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
    write.csv(filtered_sig_genes_d36e, file = file.path(joint_output_dir, "Filtered_Significant_Genes_D36E_EV_vs_MgSO4.csv"))

    # (c) Compare each active treatment against filtered D36E::EV
    for (treatment in active_treatments) {
        print(paste("  - Comparing", treatment, "vs filtered D36E::EV..."))
        res <- results(dds, contrast = c("treatment", treatment, "D36E_EV"))

        # Get normalized counts
        norm_counts <- counts(dds, normalized=TRUE)[rownames(res),]
        res$normalized_counts <- norm_counts
        
        # Remove background noise genes from active treatment results
        res_filtered <- res[!rownames(res) %in% background_genes, ]
        write.csv(as.data.frame(res_filtered), file = file.path(joint_output_dir, paste0(treatment, "_vs_Filtered_D36E_EV.csv")))
        
        # Save significant genes
        sig_genes <- subset(as.data.frame(res_filtered), padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
        write.csv(sig_genes, file = file.path(joint_output_dir, paste0("Significant_Genes_", treatment, "_vs_Filtered_D36E_EV.csv")))
        
        # Generate MA plot
        png(file.path(joint_output_dir, paste0("MA_plot_", treatment, "_vs_Filtered_D36E_EV.png")))
        plotMA(res_filtered, main = paste("MA Plot:", treatment, "vs Filtered D36E::EV"))
        dev.off()
        
        # Generate Volcano plot
        res_df <- as.data.frame(res_filtered)
        res_df$significant <- ifelse(res_df$padj < CRITICAL_VALUE & abs(res_df$log2FoldChange) > MIN_LOG2FOLD_CHANGE, "Yes", "No")
        png(file.path(working_dir, paste0("Volcano_plot_", treatment, "_vs_Filtered_D36E_EV.png")))
        print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
            geom_point() +
            theme_minimal() +
            labs(title = paste("Volcano Plot:", treatment, "vs Filtered D36E::EV"),
                x = "Log2 Fold Change",
                y = "-Log10 Adjusted P-value"))
        dev.off()
        volcano_path <- file.path(joint_output_dir, paste0("Volcano_plot_", treatment, "_", timepoint, "_vs_Filtered_D36E_EV.png"))
        ggsave(volcano_path, volcano_plot, width = 8, height = 6, dpi = 300)
    }
print("Comparisons complete.")

# # 19. Compare HopB1 vs HopAB1, excluding shared genes with D36E::EV
# print("Comparing HoB1 vs HopAB1, excluding shared genes with D36E::EV...")
# hopB1_vs_hop_AB1_output_dir <- file.path(working_dir, "/HopB1_vs_HopAB1")
# dir.create(hopB1_vs_hop_AB1_output_dir, showWarnings = FALSE)
# print("Comparing HopB1 vs HopAB1, excluding shared genes with D36E::EV...")

#     # Get significant genes expressed in D36E::EV
#     res_d36e_vs_mgso4 <- results(dds, contrast = c("treatment", "D36E_EV", "MgSO4"))
#     sig_genes_d36e <- subset(as.data.frame(res_d36e_vs_mgso4), padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
#     background_genes_d36e <- rownames(sig_genes_d36e)

#     # Perform comparison between HopB1 and HopAB1
#     res_hopb1_vs_hopab1 <- results(dds, contrast = c("treatment", "HopB1", "HopAB1"))

#     # Exclude genes shared with D36E::EV
#     filtered_res_hopb1_vs_hopab1 <- res_hopb1_vs_hopab1[!rownames(res_hopb1_vs_hopab1) %in% background_genes_d36e, ]
#     write.csv(as.data.frame(filtered_res_hopb1_vs_hopab1), file = file.path(hopB1_vs_hop_AB1_output_dir, "Filtered_HopB1_vs_HopAB1.csv"))

#     # Save significant genes for HopB1 vs HopAB1 after filtering
#     filtered_sig_genes_hopb1_vs_hopab1 <- subset(filtered_res_hopb1_vs_hopab1, padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
#     write.csv(filtered_sig_genes_hopb1_vs_hopab1, file = file.path(hopB1_vs_hop_AB1_output_dir, "Filtered_Significant_Genes_HopB1_vs_HopAB1.csv"))

#     # Generate MA plot for filtered HopB1 vs HopAB1 comparison
#     png(file.path(hopB1_vs_hop_AB1_output_dir, "MA_plot_Filtered_HopB1_vs_HopAB1.png"))
#     plotMA(filtered_res_hopb1_vs_hopab1, main = "MA Plot: Filtered HopB1 vs HopAB1")
#     dev.off()

#     # Generate Volcano plot for filtered HopB1 vs HopAB1 comparison
#     res_df <- as.data.frame(filtered_res_hopb1_vs_hopab1)
#     res_df$significant <- ifelse(res_df$padj < CRITICAL_VALUE & abs(res_df$log2FoldChange) > MIN_LOG2FOLD_CHANGE, "Yes", "No")

#     png(file.path(hopB1_vs_hop_AB1_output_dir, "Volcano_plot_Filtered_HopB1_vs_HopAB1.png"))
#     print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
#         geom_point() +
#         theme_minimal() +
#         labs(title = "Volcano Plot: Filtered HopB1 vs HopAB1",
#             x = "Log2 Fold Change",
#             y = "-Log10 Adjusted P-value"))
#     dev.off()

#     print("Comparison between HopB1 and HopAB1 complete.")


#20 TIME-RESOLVED ANALYSIS
print("Beginning time-resolved analysis with MgSO4 subtraction...")
timepoint_output_dir <- file.path(working_dir, "time_resolved_analysis")
dir.create(timepoint_output_dir, showWarnings = FALSE)
ref_timepoint <- levels(metadata$timepoint)[1] # Get reference timepoint (first level)
# First, get the background genes from MgSO4 vs D36E::EV comparison
res_d36e_vs_mgso4 <- results(dds, contrast = c("treatment", "D36E_EV", "MgSO4"))
background_genes <- rownames(subset(as.data.frame(res_d36e_vs_mgso4), padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE))
for (timepoint in levels(metadata$timepoint)) {
    for (treatment in active_treatments) {
        # Create contrast name components
        treatment_var <- paste0("treatment", treatment)
        d36e_var <- "treatmentD36E_EV"
        time_var <- paste0("timepoint", timepoint)
        
        if (timepoint == ref_timepoint) {
            # For reference timepoint, use main effect contrast
            contrast <- c("treatment", treatment, "D36E_EV")
        } else {
            # For non-reference timepoints, use interaction terms
            contrast <- list(
                c(paste0(treatment_var, ".", time_var)),
                c(paste0(d36e_var, ".", time_var))
            )
        }
        
        # Get results
        res <- results(dds, contrast = contrast)

        # Get normalized counts
        norm_counts <- counts(dds, normalized=TRUE)[rownames(res),]
        res$normalized_counts <- norm_counts
        
        # Remove background noise genes
        res_filtered <- res[!rownames(res) %in% background_genes, ]
        
        # Save filtered results
        out_file <- paste0(treatment, "_", timepoint, "_vs_Filtered_D36E_EV.csv")
        write.csv(as.data.frame(res_filtered), file.path(timepoint_output_dir, out_file))
        
        # Save significant genes
        sig_genes <- subset(as.data.frame(res_filtered), 
                          padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
        sig_file <- paste0("Significant_Genes_", treatment, "_", timepoint, "_vs_Filtered_D36E_EV.csv")
        write.csv(sig_genes, file.path(timepoint_output_dir, sig_file))

        # Generate MA plot
        ma_plot_path <- file.path(timepoint_output_dir, paste0("MA_plot_", treatment, "_", timepoint, "_vs_Filtered_D36E_EV.png"))
        png(ma_plot_path, width = 800, height = 600)
        plotMA(res_filtered, main = paste("MA Plot:", treatment, timepoint, "vs Filtered D36E_EV"))
        dev.off()
        
        # Generate Volcano plot
        volcano_data <- as.data.frame(res_filtered)
        volcano_data$significant <- ifelse(
            volcano_data$padj < CRITICAL_VALUE & abs(volcano_data$log2FoldChange) > MIN_LOG2FOLD_CHANGE,
            "Significant", "Not significant"
        )
        volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("Not significant" = "grey50", "Significant" = "red")) +
            theme_minimal(base_size = 14) +
            labs(title = paste("Volcano Plot:", treatment, timepoint, "vs Filtered D36E_EV"),
                 x = "Log2 Fold Change",
                 y = "-Log10 Adjusted P-value") +
            geom_hline(yintercept = -log10(CRITICAL_VALUE), linetype = "dashed") +
            geom_vline(xintercept = c(-MIN_LOG2FOLD_CHANGE, MIN_LOG2FOLD_CHANGE), linetype = "dashed")
        
        volcano_path <- file.path(timepoint_output_dir, paste0("Volcano_plot_", treatment, "_", timepoint, "_vs_Filtered_D36E_EV.png"))
        ggsave(volcano_path, volcano_plot, width = 8, height = 6, dpi = 300)
    }
}
print("Time-resolved analysis with MgSO4 subtraction complete.")


# 21. WITHIN-TREATMENT TIMEPOINT COMPARISONS
print("Adding within-treatment time comparisons (8h vs 1h)...")
within_treatment_dir <- file.path(working_dir, "within_treatment_time_comparisons")
dir.create(within_treatment_dir, showWarnings = FALSE)
# Verify coefficient names first
coef_names <- resultsNames(dds)
time_main <- grep("timepoint.*8h.*1h", coef_names, value = TRUE)[1]
for (treatment in active_treatments) {
    # Construct interaction term name
    interaction_term <- paste0("treatment", treatment, ".timepoint8h")
    
    if (!interaction_term %in% coef_names) {
        warning(paste("Skipping", treatment, "- interaction term not found:", interaction_term))
        next
    }
    
    # Create contrast combining main time effect + treatment-specific interaction
    contrast <- list(
        c(time_main, interaction_term),  # 8h effect components
        character(0)                     # 1h baseline
    )
    
    # Get results
    res <- results(dds, contrast = contrast)
    
    # Save results
    out_file <- paste0(treatment, "_8h_vs_1h.csv")
    write.csv(as.data.frame(res), file.path(within_treatment_dir, out_file))
    
    # Save significant genes
    sig_genes <- subset(as.data.frame(res), 
                      padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
    sig_file <- paste0("Significant_Genes_", treatment, "_8h_vs_1h.csv")
    write.csv(sig_genes, file.path(within_treatment_dir, sig_file))
    
    # Generate MA plot
    png(file.path(within_treatment_dir, paste0("MA_plot_", treatment, "_8h_vs_1h.png")))
    plotMA(res, main = paste("MA Plot:", treatment, "8h vs 1h"))
    dev.off()
    
    # Generate Volcano plot
    volcano_df <- as.data.frame(res)
    volcano_df$significant <- ifelse(
        volcano_df$padj < CRITICAL_VALUE & abs(volcano_df$log2FoldChange) > MIN_LOG2FOLD_CHANGE,
        "Yes", "No"
    )
    png(file.path(within_treatment_dir, paste0("Volcano_plot_", treatment, "_8h_vs_1h.png")))
    print(ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
        geom_point() +
        scale_color_manual(values = c("No" = "grey60", "Yes" = "red")) +
        theme_minimal() +
        labs(title = paste(treatment, "8h vs 1h"),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"))
    dev.off()
    volcano_path <- file.path(within_treatment_dir, paste0("Volcano_plot_", treatment, "_", timepoint, "_8h_vs_1h.png"))
    ggsave(volcano_path, volcano_plot, width = 8, height = 6, dpi = 300)
}
print("Within-treatment time comparisons complete.")
