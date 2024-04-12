library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')


#'
#'timepoint from sample function from last assignment
#'used as helper function
#'
timepoint_from_sample <- function(x) {
  find <- str_locate(x, "_") 
  timepoint <- str_sub(x, 1, find -1)
  timepoint <- unique(timepoint)
  return(factor(timepoint))
}

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(csv_path, metafile, selected_times) {
  counts_data <- read.csv(csv_path, sep = "\t")
  meta_data <- read.csv(metafile)

  rowData <- counts_data["gene"]
  counts_data <- as.data.frame(dplyr::select(counts_data, -gene))
  rownames(counts_data) <- rowData$gene
  
  colData <- DataFrame(
    samplename = colnames(counts_data),
    timepoint = sapply(colnames(counts_data), timepoint_from_sample)
  )
  
  colData <- colData[colData$timepoint %in% selected_times, ]
  
  selected_samples <- colData$samplename
  counts_data <- counts_data[, selected_samples]

  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts_data)),
    colData = colData,
  )
  
  metadata(se) <- list(model = "model")
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  unique_genes <- rownames(se)[!duplicated(rownames(se))]
  se_unique <- se[unique_genes, ]
  
  dds <- DESeqDataSet(se_unique, design = design)
  dds <- DESeq(dds)
  results <- results(dds, contrast = c("timepoint", "vAd", "vP0"))
  results_df <- as.data.frame(
    results
  )
  result_list <- list(results_df = results_df, dds= dds)
  return(result_list)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  
  deseq2_results <- deseq2_res$results_df
  
  labeled <- as_tibble(rownames = "genes", deseq2_results) %>%
    mutate(volc_plot_status = case_when(
      deseq2_results$padj < padj_threshold & deseq2_results$log2FoldChange > 0 ~ "UP",
      deseq2_results$padj < padj_threshold & deseq2_results$log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    )) 
  
  labeled_with_genes <- labeled %>% mutate(genes = rownames(deseq2_results))
  
  labeled_final <- labeled_with_genes %>% 
    dplyr::select(genes, volc_plot_status, everything())
  
  return(labeled_final)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_hist <- ggplot(labeled_results, aes(x = pvalue)) + 
    geom_histogram(binwidth = 0.02, fill = "lightblue", color = "black") + 
    labs(
      x = "pvalue" ,
      y = "counts" ,
      title = "Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)"
    ) + 
    theme_light()
  return(pval_hist)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_genes <- labeled_results %>%
    filter(padj < padj_threshold)
  
  log2fc_hist <- ggplot(sig_genes, aes(x = log2FoldChange)) + 
    geom_histogram(binwidth = 0.2, fill = "lightblue", color = "black") + 
    labs(
      x = "log2FoldChange",
      y = "count",
      title = "Histogram of Log2FoldChange for DE genes (vP0 vs. vAd)"
    ) +
    theme_light()
  return(log2fc_hist)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- labeled_results %>%
    arrange(padj) %>%
    head(num_genes)
  
  top_gene_names <- top_genes$genes
  dds_subset <- dds_obj[top_gene_names, ]
  
  normalized_counts <- counts(dds_subset, normalized = TRUE)
  
  dds_merged <- as.data.frame(normalized_counts) %>%
    rownames_to_column(var = "genes") %>%
    pivot_longer(-genes, names_to = "sample", values_to = "counts")
  
  scatter_plot <- ggplot(dds_merged, aes(x = genes, y = log10(counts), color = sample)) +
    geom_point(position = position_jitter(width = 0.2)) +
    scale_y_continuous(trans = "log10") +
    labs(
      x = "Genes",
      y = "Normalized Counts",
      title = "Plot of log10(normalized counts) from top 10 DE genes"
    ) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(scatter_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  volc_plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point() +
    labs(x = "log2FoldChange",
         y = "-log10(padj)",
         title = "Volcano plot of DESeq2 differential expression results vP0 vs. vAd)"
         ) + 
    scale_color_manual(values = c("UP" = "cornflowerblue", "DOWN" = "coral", "NS" = "forestgreen")) +
    theme_light() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  return(volc_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene_map <- read.table(id2gene_path, header = FALSE)
  colnames(id2gene_map) <- c("EnsemblID", "MGI_Symbol")
  
  merged_data <- merge(labeled_results, id2gene_map, by.x = "genes", by.y = "EnsemblID", all.x = TRUE)
  filtered_data <- filter(merged_data, !is.na(log2FoldChange))
  sorted_data <- filtered_data[order(-filtered_data$log2FoldChange), ]
  
  ranked <- setNames(sorted_data$log2FoldChange, sorted_data$MGI_Symbol)
  return(ranked)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  rnk_list <- rnk_list[!is.na(rnk_list)]
  rnk_list <- rnk_list[is.finite(rnk_list)]

  pathways <- fgsea::gmtPathways(gmt_file_path)
  fgsea_res <- fgsea(pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )
  return(tibble(fgsea_res))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
   top_positive_nes <- fgsea_results %>%
     filter(padj < .25 & NES > 0) %>%
     slice_max(NES, n = num_paths)
  
   top_negative_nes <- fgsea_results %>%
     filter(padj < 0.25 & NES < 0) %>%
     slice_min(NES, n = num_paths)

   combined <- bind_rows(top_positive_nes, top_negative_nes)
   
   combined <- combined %>%
     mutate(pathway = forcats::fct_reorder(pathway, NES),
            group = ifelse(NES > 0, "top_positive", "top_negative"))
  
  TP_plot <- combined %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
    ggplot() +
    geom_bar(aes(x = pathway, y = NES, fill = group), stat = "identity") +
    scale_fill_manual(values = c("top_positive" = "red", "top_negative" = "blue")) +
    theme_light() +
    ggtitle("fgsea results for C2 conical Pathways") +
    ylab("normalized enrichment score (NES)") +
    xlab("") +
    coord_flip()
  
  return(TP_plot)
}