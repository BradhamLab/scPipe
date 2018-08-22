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


count_data <- read.csv('output/matrix/filtered_count_matrix.csv',
                         row.names=1, check.names=FALSE)
metadata <- read.csv('output/metadata/filtered_metadata.csv',
                     row.names=1, check.names=FALSE)
#' Title
#'
#' @param count_matrix
#' @param metadata
#' @param min_cluster_size
#' @param size_step
#'
#' @return
#' @export
#'
#' @examples
normalize_data <- function(count_matrix, metadata, cluster=FALSE,
                           min_cluster_size=5, size_step=5) {
  batch_sizes <- sapply(unique(metadata$batch), function(x) {
    return(sum(metadata$batch == x))
  })
  counts = as.matrix(count_matrix)
  log_counts <- log_transform_counts(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(
           assays=list('counts'=counts, 'logcounts'=log_counts),
           colData=metadata)
  clusters <- NULL
  max_size <- size_step*floor(max(batch_sizes)/size_step)
  if (cluster) {
    clusters <- scran::quickCluster(sce, min.size=min_cluster_size,
                                    max.size=max_size, assay.type='logcounts')
    max_size = size_step*floor(min(table(clusters))/size_step)
  }
  sce <- scran::computeSumFactors(sce, size=seq(20, max_size, size_step),
                                  cluster=clusters, assay.type='counts')
  norm <- scater::normalizeSCE(sce, exprs_values='counts',
                               return_norm_as_exprs=True)
  return(norm)
}

#' Title
#'
#' @param count_matrix
#' @param zero_sub
#'
#' @return
#' @export
#'
#' @examples
log_transform_counts <- function(count_matrix, zero_sub=10^(-6)) {
  log_counts = log10(as.matrix(count_matrix))
  log_counts[log_counts == -Inf] = log10(10^(-6))
  return(log_counts)
}


#' Title
#'
#' @param sce
#'
#' @return
#' @export
#'
#' @examples
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


#' create_boxplot
#'
#' Create grouped violin/boxplots for a specified feature.
#'
#' @param dataframe (data.frame): data.frame containing data to plot.
#' @param y_column (string): column name to plot along the y-axis.
#' @param group_column (string): column name for grouping variable.
#' @param facet (string, optional): second grouping variable. Will create
#'     distinct boxplots for each unique value within the column. Boxplots are
#'     stacked in rows. Default is '', and no facet wrapping is applied.
#'
#' @return (gg.ggplot): ggplot object of boxplots
#' @export
#'
#' @examples
#'
#' # plot 'Pmi' accross disease states
#' create_boxplots(combined_data, 'Pmi', 'Disease.state')
#'
#' # plot 'Pmi' accross disease states and tissue types
#' create_boxplots(combined_data, 'Pmi', 'Disease.state', 'Tissue')
create_boxplot <- function(dataframe, y_column, group_column, facet='') {
  boxplots <- ggplot(dataframe, aes_string(x=group_column, y=y_column)) +
    geom_violin(aes_string(fill=group_column), trim=FALSE) +
    geom_boxplot(width=0.1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  if (facet != '' && facet %in% colnames(dataframe)) {
    boxplots <- boxplots + facet_grid(reformulate('.', facet))
  }
  return(boxplots)
}


visualize_normalization <- function(counts, norm) {
  # find most extreme cells in count data
  cell_counts <- sapply(colnames(counts), function(x) {
    return(sum(counts[ , x]))})
  ordered_cells <- order(cell_counts, decreasing=TRUE)
  n <- length(ordered_cells)
  top_5 <- names(cell_counts[ordered_cells[1:5]])
  bottom_5 <- names(cell_counts[ordered_cells[n:(n-4)]])

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
  norm_df <- as.data.frame(t(exprs(norm)[ , c(top_5, bottom_5)]))
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
  p <- ggplot(combined, aes(x=Cell, y=Expression)) +
              geom_violin(aes(fill=Count.Level), trim=FALSE) +
              geom_boxplot(width=0.1) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle=45, hjust=1)) +
              facet_grid(~Normalized)
}

visualize_batch_effect <- function(norm, norm_batch) {
  batch_umap <- as.data.frame(umap(t(logcounts(norm)))$layout)
  colnames(batch_umap) <- c('UMAP1', 'UMAP2')
  no_batch_umap <- as.data.frame(umap(t(logcounts(norm_batch)))$layout)
  colnames(no_batch_umap) <- c('UMAP1', 'UMAP2')

  batch_umap[row.names(colData(norm)), 'batch'] <- colData(norm)$batch
  batch_umap$Batch.Normalized <- 'Before'
  no_batch_umap[row.names(colData(norm_batch)), 'batch'] <- colData(norm_batch)$batch
  no_batch_umap$Batch.Normalized <- 'After'
  combined <- rbind(batch_umap, no_batch_umap)
  p <- ggplot(combined, aes(x=UMAP1, y=UMAP2, col=batch)) +
         geom_point() +
         ggtitle('Before and After Batch Correction') +
         theme(legend.title=element_text('Batch')) +
         facet_wrap(~Batch.Normalized)

  return(p)
}

main <- function(count_data, metadata) {
  norm_data <- normalize_data(count_data, metadata, cluster=FALSE)
  batch_norm <- remove_batch_effects(norm_data)
  visualize_batch_effect(norm_data, batch_norm)
}


