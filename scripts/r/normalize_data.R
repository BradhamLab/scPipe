#' Script to normalize and batch-correct single-cell RNAseq data.
#'
#' @author: Dakota Hawkins
#' @date: August 20, 2018


suppressMessages(require(scran))
suppressMessages(require(scater))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(monocle))
suppressMessages(require(scater))
suppressMessages(require(umap))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(cowplot))

#' Normalize Data
#'
#' Normalize count data via `scran` normalization using cell
#' pooling.
#'
#' @param count_matrix (data.frame): gene x sample matrix of
#'   read counts.
#' @param metadata (data.frame): data.frame containing cell
#'   meta data including a `batch` column for cell batches.
#' @export
batch_normalize_data <- function(count_matrix, metadata) {

  # create SingleCellExperiment object for normalization
  counts <- as.matrix(count_matrix)
  log_counts <- log_transform_counts(counts)
  sce <- SingleCellExperiment::SingleCellExperiment(
           assays=list('counts'=counts, 'logcounts'=log_counts),
           colData=metadata)

  # compute batch-specific sum factors
  batches <- unique(metadata$batch)
  samples_per_batch <- sapply(batches, function(x) {
    sum(metadata$batch == x)
  })
  by_batch <- lapply(batches[-samples_per_batch], function(x) {
    return(scran::computeSumFactors(filter(sce, batch == x)))
    })
  # set multiBatchNorm arguments
  by_batch$norm.args <- list('return_log' = TRUE)
  batch_norm <- do.call(scran::multiBatchNorm, by_batch)

  # combine batch normalized sce's
  norm <- do.call(cbind, batch_norm)
  return(norm)
}


#' Log Transform Counts Data
#'
#' Log2 transform count data.
#'
#' Log2 transform count data while providing an offset to avoid negative and
#' undefined values.
#'
#' @param count_matrix (data.frame): gene x sample read count data.
#' @param offset (float): value to displace values by to ensure non-negative
#'   log values. Default is 1.
#'
#' @return (data.frame): log-transformed read counts.
#' @export
log_transform_counts <- function(count_matrix, offset=1) {
  return(log2(as.matrix(count_matrix) + offset))
}

#' Estimate absolute abundance of transcripts
#'
#' Estimates absolute abundance of transcripts from TPM data using the Census
#' Method outlined in Monocle.
#'
#' @param tpm_matrix (data.frame): g x sample tpm data
#'
#' @return (data.matrix): matrix of estimate counts.
#'
#' @references
#'   Qiu, X., Hill, A., Packer, J., Lin, D., Ma, Y.-A., & Trapnell, C. (2017).
#'     Single-cell mRNA quantification and differential analysis with Census.
#'     Nature Methods, 14(3), 309–315. https://doi.org/10.1038/nmeth.4150
absolute_abundance <- function(tpm_matrix) {
  tpm_cds <- suppressWarnings(monocle::newCellDataSet(as.matrix(tpm_matrix)))
  abs_matrix <- monocle::relative2abs(tpm_cds, method='tpm_fraction')
  return(abs_matrix)
}

#' Visualize Normalization
#'
#' Plot expression distributions for the five cells with the highest and
#' lowest library sizes in original count data.
#'
#' @param counts (data.frame): gene x sample read counts.
#' @param norm (SingleCellExperiment): SingleCellExperiment with
#'   normalized count data in the `exprs` slot.
#'
#' @return (ggplot): violin + boxplots showing differences in cells
#'  pre and post normalization.
#' @export
visualize_normalization <- function(counts, norm) {
  # find most extreme cells in count data
  cell_counts <- sapply(colnames(counts), function(x) {
    return(sum(counts[ , x]))})
  ordered_cells <- order(cell_counts, decreasing=TRUE)
  n <- length(ordered_cells)
  top_5 <- names(cell_counts[ordered_cells[1:5]])
  bottom_5 <- names(cell_counts[ordered_cells[n:(n-4)]])

  cell_order <- c(sort(top_5), sort(bottom_5))

  # log transform counts in extreme cells
  log_count_df <- as.data.frame(
                    log_transform_counts(t(counts[ , c(top_5, bottom_5)])))
  zero_measures <- log_count_df == 0
  log_count_df[zero_measures] <- NA
  log_count_df[top_5, 'Count.Level'] <- 'High'
  log_count_df[bottom_5, 'Count.Level'] <- 'Low'
  log_count_df$Normalized <- 'Before'
  log_count_df$Cell <- row.names(log_count_df)

  # make long format
  id_vars <- c('Cell', 'Count.Level', 'Normalized')
  melted_count <- melt(log_count_df, id.vars=id_vars, variable.name='Gene',
                       value.name='Expression')

  # get extreme log values in normalized data
  norm_df <- as.data.frame(t(logcounts(norm)[ , c(top_5, bottom_5)]))
  norm_df[zero_measures] <- NA
  norm_df[top_5, 'Count.Level'] <- 'High'
  norm_df[bottom_5, 'Count.Level'] <- 'Low'
  norm_df$Normalized <- 'After'
  norm_df$Cell <- row.names(norm_df)

  # make long format
  melted_norm <- melt(norm_df, id.vars=id_vars, variable.name='Gene',
                       value.name='Expression')

  # combine data for plotting
  combined <- rbind(melted_count, melted_norm)
  combined$Cell <- factor(combined$Cell, levels=cell_order)
  combined$Normalized <- factor(combined$Normalized, levels=c('Before',
                                                              'After'))
  p <- ggplot(combined, aes(x=Cell, y=Expression)) +
              geom_violin(aes(fill=Count.Level), trim=FALSE) +
              geom_boxplot(width=0.1) +
              theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
              facet_grid(~Normalized)
  return(p)
}


#' Visualize Batch Effect Removal
#'
#' Plot cells on UMAP dimensions, coloring by batch, to visualize
#' batch dispersion pre and post batch effect removal.
#'
#' @param norm (SingleCellExperiment): normalized data before batch
#'  effect removal
#' @param norm_batch (SingleCellExperiment): normalized data after
#'  batch removal
#'
#' @return (ggplot): scatter plot visualizing batch dispersion
#'  pre and post batch effect removal
#' @export
visualize_batch_effect <- function(counts, norm_batch) {
  batch_umap <- as.data.frame(umap(t(log_transform_counts(counts)))$layout)
  colnames(batch_umap) <- c('UMAP1', 'UMAP2')
  no_batch_umap <- as.data.frame(umap(t(logcounts(norm_batch)))$layout)
  colnames(no_batch_umap) <- c('UMAP1', 'UMAP2')

  batch_umap[row.names(colData(norm_batch)), 'Batch'] <- colData(norm_batch)$batch
  batch_umap$Batch.Normalized <- 'Before'
  no_batch_umap[row.names(colData(norm_batch)), 'Batch'] <- colData(norm_batch)$batch
  no_batch_umap$Batch.Normalized <- 'After'
  combined <- rbind(batch_umap, no_batch_umap)
  combined$Batch.Normalized <- factor(combined$Batch.Normalized,
                                      c('Before', 'After'))
  p <- ggplot(combined, aes(x=UMAP1, y=UMAP2, col=Batch)) +
         geom_point() +
         ggtitle('Before and After Batch Correction') +
         facet_wrap(~Batch.Normalized)

  return(p)
}

#' Preprocess Read Count Expression Data
#'
#' @param count_data (data.frame): gene x sample read counts.
#' @param tpm_data (data.frame): gene x sample tpm matrix.
#' @param metadata (data.frame): metadata containing information for cells.
#'   including a `batch` column.
#'
#' @return (list): named list containing between-sample and batch normalized
#'   data ('data'), between-sample visualization ('norm_viz'), and batch
#'   correction visualization ('batch_viz').
#'
#' @export
main <- function(count_data, tpm_data, metadata) {
  out_data <- lapply(list(count_data, tpm_data), function(data) {
    normalized <- batch_normalize_data(data, metadata)
    norm_viz <- visualize_normalization(data, normalized)
    batch_viz <- visualize_batch_effect(data, normalized)
    return(list('data'=normalized, 'norm_viz'=norm_viz, 'batch_viz'=batch_viz))
  })
  names(out_data) <- c('count', 'tpm')
  return(out_data)
}

if (exists('snakemake')) {
  # read in data
  count_data <- read.csv(snakemake@input[['cmat']],
                         row.names=1, check.names=FALSE)
  tpm_data <- read.csv(snakemake@input[['tpm']],
                       row.names=1, check.names=FALSE)
  metadata <- read.csv(snakemake@input[['meta']],
                       row.names=1, check.names=FALSE)

  # check if plot dir exists, create if no
  if (!dir.exists(snakemake@params[['plot_dir']])) {
    dir.create(snakemake@params[['plot_dir']])
  }

  # get processed data and visualizations
  out_data <- main(count_data, tpm_data, metadata)

  # write normalized log-count data to 'final' directory
  write.csv(logcounts(out_data$count$data), file=snakemake@output[['cmat']],
            quote=FALSE)

  # write normalized tpm data to 'final' directory
  write.csv(logcounts(out_data$tpm$data), file=snakemake@output[['tpm']],
            quote=FALSE)

  # write metadata to 'final' directory
  write.csv(metadata, file=snakemake@output[['meta']],
            quote=FALSE)

  # write normalization plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'count_normalization.png'),
         plot=out_data$count$norm_viz,
         device='png')

  # write batch-effect removal plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'count_batch_removal.png'),
         plot=out_data$count$batch_viz,
         device='png')

  # write normalization plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'tpm_normalization.png'),
         plot=out_data$tpm$norm_viz,
         device='png')

  # write batch-effect removal plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'tpm_batch_removal.png'),
         plot=out_data$tpm$batch_viz,
         device='png')
}
