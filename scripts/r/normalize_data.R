#' Script to normalize and batch-correct single-cell RNAseq data.
#'
#' @author: Dakota Hawkins
#' @date: August 20, 2018


require(scran)
require(scater)
require(SingleCellExperiment)
require(scater)
require(umap)
require(ggplot2)
require(reshape2)


#' Normalize Data
#'
#' Normalize count data via `scran` normalization using cell
#' pooling.
#'
#' @param count_matrix (data.frame): gene x sample matrix of
#'   read counts.
#' @param metadata (data.frame): data.frame containing cell
#'   meta data including a `batch` column for cell batches.
#' @param cluster (boolean): whether to perform clustering
#'   prior to cell pooling.
#' @param min_cluster_size (int): minimum number of cells
#'   per cluster.
#' @param step_size (int): value to increment step by when
#'   setting the size sequence in `scran::ComputeSumFactors`.
#'
#' @return (SingleCellExperiment): SingleCellExperiment with
#'   normalized expression values in the `exprs` slot.
#' @export
normalize_data <- function(count_matrix, metadata, cluster=FALSE,
                           min_cluster_size=5, step_size=5) {
  batch_sizes <- sapply(unique(metadata$batch), function(x) {
    return(sum(metadata$batch == x))
  })
  counts = as.matrix(count_matrix)
  log_counts <- log_transform_counts(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(
           assays=list('counts'=counts, 'logcounts'=log_counts),
           colData=metadata)
  clusters <- NULL
  max_size <- step_size * floor(ncol(count_matrix)/2/step_size)
  if (cluster) {
    clusters <- scran::quickCluster(sce, min.size=min_cluster_size,
                                    max.size=max_size, assay.type='logcounts')
    max_size = step_size*floor(min(table(clusters))/step_size)
  }
  sce <- scran::computeSumFactors(sce, size=seq(20, max_size, step_size),
                                  cluster=clusters, assay.type='counts')
  norm <- scater::normalizeSCE(sce, exprs_values='counts', return_log=FALSE,
                               return_norm_as_exprs=True)
  logcounts(norm) <- log_transform_counts(normcounts(norm))
  return(norm)
}

#' Log Transform Counts Data
#'
#' Log10 transform count data. Sets zeros to `zero_sub` to avoid negative
#' infinities.
#'
#' @param count_matrix (data.frame): gene x sample read count data.
#' @param zero_sub (float): replacement value for zeros.
#'
#' @return (data.frame): log-transformed read counts.
#' @export
log_transform_counts <- function(count_matrix, zero_sub=10^(-6)) {
  log_counts = log10(as.matrix(count_matrix))
  log_counts[log_counts == -Inf] = log10(10^(-6))
  return(log_counts)
}


#' Remove Batch Effects
#'
#' Remove batch effects from expression data.
#'
#' @param sce (SingleCellExperiment): SingleCellExperiment object
#'   with normalized count data.
#'
#' @return (SingleCellExperiment): batch corrected data.
#' @export
remove_batch_effects <- function(sce, n_mnn=20) {
  # subset log10 normalized expression matrix by batch
  batches <- levels(colData(sce)$batch)
  expr_batches <- lapply(batches, function(x) {
    cells <- row.names(subset(colData(sce), batch == x))
    expr_data <- logcounts(sce)[ , cells]
    return(expr_data)
  })
  # make largest batch reference batch
  cells_per_batch <- sapply(expr_batches, ncol)
  expr_batches <- expr_batches[order(cells_per_batch, decreasing=TRUE)]

  # set mnnCorrect arguments
  expr_batches$k <- n_mnn
  expr_batches$cos.norm.in <- TRUE
  expr_batches$cos.norm.out <- FALSE  # can change this depending on downstream
  corrected <- do.call(scran::mnnCorrect, expr_batches)

  # combine corrected expression matrices
  combined <- do.call(cbind, corrected$corrected)
  combined <- combined[ , colnames(logcounts(sce))] # maintain original order

  # create batch-corrected sce
  out <- SingleCellExperiment::SingleCellExperiment(
           assays=list('logcounts'=combined),
           colData=colData(sce))
  return(out)
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
  zero_measures <- log_count_df == -6
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
visualize_batch_effect <- function(norm, norm_batch) {
  batch_umap <- as.data.frame(umap(t(logcounts(norm)))$layout)
  colnames(batch_umap) <- c('UMAP1', 'UMAP2')
  no_batch_umap <- as.data.frame(umap(t(logcounts(norm_batch)))$layout)
  colnames(no_batch_umap) <- c('UMAP1', 'UMAP2')

  batch_umap[row.names(colData(norm)), 'Batch'] <- colData(norm)$batch
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
#' @param count_data (data.frame): gene x sample read counts
#' @param metadata (data.frame): metadata containing information for cells.
#'   including a `batch` column.
#'
#' @return (list): named list containing between-sample and batch normalized
#'   data ('data'), between-sample visualization ('norm_viz'), and batch
#'   correction visualization ('batch_viz').
#'
#' @export
main <- function(count_data, metadata) {
  norm_data <- normalize_data(count_data, metadata, cluster=FALSE)
  norm_viz <- visualize_normalization(count_data, norm_data)
  batch_norm <- remove_batch_effects(norm_data)
  batch_viz <- visualize_batch_effect(norm_data, batch_norm)
  out_list <- list('data'=batch_norm, 'norm_viz'=norm_viz,
                   'batch_viz'=batch_viz)
  return(out_list)
}

if (exists('snakemake')) {
  # read in data
  count_data <- read.csv(snakemake@input[['cmat']],
                         row.names=1, check.names=FALSE)
  metadata <- read.csv(snakemake@input[['meta']],
                       row.names=1, check.names=FALSE)

  # check if plot dir exists, create if no
  if (!dir.exists(snakemake@output[['plot_dir']])) {
    dir.create(snakemake@output[['plot_dir']])
  }

  # get processed data and visualizations
  out_data <- main(count_data, metadata)

  # write normalized log data to 'final' directory
  write.csv(logcounts(out_data$data), file=snakemake@output[['cmat']],
            quote=FALSE)

  # write metadata to 'final' directory
  write.csv(colData(out_data$data), file=snakemake@output[['meta']],
            quote=FALSE)

  # write normalization plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'normalization.png'), plot=out_data$norm_viz,
         device='png', width=10, height=8, dpi=600, units='in')

  # write batch-effect removal plot
  ggsave(filename=file.path(snakemake@params[['plot_dir']],
                            'batch_removal.png'), plot=out_data$batch_viz,
         device='png', width=10, height=8, dpi=600, units='in')
}
