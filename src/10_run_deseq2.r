library(DESeq2)
library(ggplot2)
library(org.At.tair.db)
library(AnnotationDbi)

# Constants
CRITICAL_VALUE <- 0.01
MIN_LOG2FOLD_CHANGE <- 2

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
save(dds, file = file.path(working_dir, "deseq_object.RData"))

# 18. Timepoint comparisons
print("Performing corrected timepoint comparisons...")
# Get reference timepoint and active treatments
ref_timepoint <- levels(metadata$timepoint)[1]  # Should be "1h"
timepoints <- levels(metadata$timepoint)
active_treatments <- levels(metadata$treatment)
active_treatments <- active_treatments[!active_treatments %in% c("MgSO4", "D36E_EV")]
for (timepoint in timepoints) {
    timepoint_dir <- file.path(working_dir, paste0("timepoint_", timepoint))
    dir.create(timepoint_dir, showWarnings = FALSE)
    
    # (A) D36E_EV vs MgSO4 at this timepoint
    print(paste("Comparing D36E_EV vs MgSO4 at", timepoint))
    if(timepoint == ref_timepoint) {
        res_d36e_vs_mgso4 <- results(dds, name = "treatment_D36E_EV_vs_MgSO4")
    } else {
        interaction_term <- paste0("treatmentD36E_EV.timepoint", timepoint)
        res_d36e_vs_mgso4 <- results(dds, contrast = list(
            c("treatment_D36E_EV_vs_MgSO4", interaction_term)
        ))
    }
    # Save background genes for this timepoint
    sig_genes_d36e <- subset(as.data.frame(res_d36e_vs_mgso4), 
                           padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE)
    background_genes <- rownames(sig_genes_d36e)
    

    # (B) Active treatments vs D36E_EV at this timepoint
    for (treatment in active_treatments) {
        print(paste("Comparing", treatment, "vs D36E_EV at", timepoint))
        
        # Build contrast
        if(timepoint == ref_timepoint) {
            res <- results(dds, contrast = list(
                c(paste0("treatment_", treatment, "_vs_MgSO4")),
                c("treatment_D36E_EV_vs_MgSO4")
            ))
        } else {
            main_effect <- paste0("treatment_", treatment, "_vs_MgSO4")
            interaction_term <- paste0("treatment", treatment, ".timepoint", timepoint)
            res <- results(dds, contrast = list(
                c(main_effect, interaction_term),
                c("treatment_D36E_EV_vs_MgSO4", paste0("treatmentD36E_EV.timepoint", timepoint))
            ))
        }
        
        # Filter and annotate
        res_filtered <- res[!rownames(res) %in% background_genes, ]
        sig_genes <- subset(as.data.frame(res_filtered),
                          padj < CRITICAL_VALUE & 
                          abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE &
                          !is.na(padj))
        
        # Add gene symbols
        sig_genes$symbol <- mapIds(org.At.tair.db,
                                 keys = rownames(sig_genes),
                                 column = "SYMBOL",
                                 keytype = "TAIR",
                                 multiVals = "first")
        
        # Save results
        write.csv(as.data.frame(res_filtered),
                file.path(timepoint_dir, paste0(treatment, "_vs_D36E_EV_", timepoint, ".csv")))
        
        write.csv(sig_genes,
                file.path(timepoint_dir, paste0("SigGenes_", treatment, "_vs_D36E_EV_", timepoint, ".csv")))
        
        # Generate MA/Volcano plots
        ma_plot_path <- file.path(timepoint_dir, paste0("MA_", treatment, "_", timepoint, ".png"))
        png(ma_plot_path)
        plotMA(res_filtered, main = paste(treatment, "vs D36E_EV at", timepoint))
        dev.off()
    }
}
print("Timepoint comparisons complete.")





# 19. Within-treatment timepoint comparisons
print("Adding within-treatment time comparisons (8h vs 1h)...")
within_treatment_dir <- file.path(working_dir, "within_treatment_time_comparisons")
dir.create(within_treatment_dir, showWarnings = FALSE)
# Verify coefficient names first
coef_names <- resultsNames(dds)
if (!"timepoint_8h_vs_1h" %in% coef_names) {
    stop("Main time effect (timepoint_8h_vs_1h) not found in resultsNames")
}
for (treatment in active_treatments) {
    # Construct interaction term name
    interaction_term <- paste0("treatment", treatment, ".timepoint8h")
    
    # Check if required coefficients exist
    if (!interaction_term %in% coef_names) {
        warning(paste("Skipping", treatment, "- interaction term not found:", interaction_term))
        next
    }
    
    # Create contrast combining main time effect + treatment-specific interaction
    contrast <- list(
        c("timepoint_8h_vs_1h", interaction_term),
        character(0)  # Baseline (1h) is implicit
    )
    
    # Get results
    res <- results(dds, contrast = contrast)
    
    # Filter significant genes
    sig_genes <- subset(as.data.frame(res), 
                      padj < CRITICAL_VALUE & abs(log2FoldChange) > MIN_LOG2FOLD_CHANGE & !is.na(padj))
    
    if (nrow(sig_genes) == 0) {
        warning(paste("No significant genes found for", treatment, "at 8h vs 1h"))
        next
    }
    
    # Save results
    out_file <- paste0(treatment, "_8h_vs_1h.csv")
    write.csv(as.data.frame(res), file.path(within_treatment_dir, out_file))
    
    # Save significant genes
    sig_file <- paste0("Significant_Genes_", treatment, "_8h_vs_1h.csv")
    write.csv(sig_genes, file.path(within_treatment_dir, sig_file))
    
    # Generate Volcano plot
    volcano_df <- as.data.frame(res)
    volcano_df$significant <- ifelse(
        volcano_df$padj < CRITICAL_VALUE & abs(volcano_df$log2FoldChange) > MIN_LOG2FOLD_CHANGE,
        "Yes", "No"
    )
    
    volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
        geom_point() +
        scale_color_manual(values = c("No" = "grey60", "Yes" = "red")) +
        theme_minimal() +
        labs(title = paste(treatment, "8h vs 1h"),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value") +
        geom_hline(yintercept = -log10(CRITICAL_VALUE), linetype = "dashed") +
        geom_vline(xintercept = c(-MIN_LOG2FOLD_CHANGE, MIN_LOG2FOLD_CHANGE), linetype = "dashed")
    
    volcano_path <- file.path(within_treatment_dir, paste0("Volcano_plot_", treatment, "_8h_vs_1h.png"))
    ggsave(volcano_path, volcano_plot, width = 8, height = 6, dpi = 300)
}
print("Within-treatment time comparisons complete.")
