#' Script to normalize and batch-correct single-cell RNAseq data.
#'
#' @author: Dakota Hawkins
#' @date: August 20, 2018


require(scran)
require(scater)
require(ggplot2)

count_matrix <- read.csv('output/matrix/filtered_count_matrix.csv',
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
normalize_data <- function(count_matrix, metadata, min_cluster_size=50,
                           size_step=5) {
  sce <- SingleCellExperiment(assays=list('counts'=as.matrix(count_matrix)),
                              colData=metadata)

  clusters <- scran::quickCluster(sce, min.size=min_cluster_size)
  max_size = size_step*floor(min(table(clusters))/size_step)
  sce <- computeSumFactors(sce, size=seq(size_step, max_size, size_step),
                           cluster=clusters)
  norm <- scater::normalize(sce)
  return(norm)
}

remove_batch_effects <- function(sce) {

  return(NULL)
}

visualize_batch_effect <- function() {

  return(NULL)
}

