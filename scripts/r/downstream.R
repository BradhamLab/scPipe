library(scran)
library(scater)
library(umap)
library(ggplot2)
library(cowplot)
library(monocle)
library(GSEABase)
library(GSVA)
library(inflection)

#             Lv-Ets1/2     Lv-Msp130   Lv-Pdgfr/vegfr   Lv-Sm30b
pmc_spus <-c('SPU_002874', 'SPU_002088', 'SPU_000310', 'SPU_000826',
             #  Lv-Sm50     Lv-erg         Lv-Fox0_1    LV-pks2
             'SPU_018811',  'SPU_018483', 'SPU_028698', 'SPU_028395',
             # tbx2/3       frp (FrED)
             'SPU_023386', 'SPU_022598')

subset_genes <- c('SPU_022598', 'SPU_004767', 'SPU_000826', 'SPU_003102',
                  'SPU_028395', 'SPU_023386')


tpm_data <- read.csv('/home/dakota/PythonModules/scPipe/output/final/normalized_tpm_matrix.csv',
                     row.names=1, check.names=FALSE)
count_data <- read.csv('/home/dakota/PythonModules/scPipe/output/final/normalized_log_matrix.csv',
                       row.names=1, check.names=FALSE)
cell_data <- read.csv('/home/dakota/PythonModules/scPipe/output/metadata/filtered_metadata.csv',
                      row.names=1, check.names=FALSE)
cell_data$cell.id <- row.names(cell_data)
gene_data <- read.csv('/home/dakota/SequenceData/evm_annotations.csv',
                      row.names=1, check.names=FALSE)

bad_controls <- c("ASW_C07_2016-07-29", "ASW_C01_2016-07-29",
                  "ASW_C09_2016-07-29", "ASW_A02_2016-07-29",
                  "ASW_G06_2016-07-29", "ASW_E05_2016-07-29",
                  "ASW_A12_2016-07-29", "ASW_G03_2016-07-29",
                  "ASW_E04_2016-07-29", "ASW_E11_2016-07-29",
                  "ASW_C02_2016-07-29")



#' Create a `GSEABase` gene set object
#'
#' @param gene_vector (vector, char): vector of gene ids that form the gene
#'   set.
#' @param name (char): name of gene set.
#'
#' @return A `GSEABase` gene set object.
create_gene_set <- function(gene_vector, name) {
  gene_set <- GSEABase::GeneSet(gene_vector, setName = name)
  return(gene_set)
}

#' Run enrichment analysis
#'
#' @param e_set (`BioBase:ExpressionSet`): `ExpressionSet` containing
#'   gene-expression profiles.
#'
#' @param gene_collect (`GSEABase:GeneSetCollection`): `GeneSetCollection`
#'   containing marker genes for phenotype of interest
#'
#' @param method (char): method for enrichment analysis. See `GSVA`
#'   documentaion for more information.
#'
#' @return An (m x n) data matrix where m is the number of phenotypes of
#'   interest, and n is the number of samples in `e_set`.
enrichment <- function(e_set, gene_collection, method = 'ssgsea') {
  enrichment <- GSVA::gsva(exprs(e_set), gene_collection, method = method)
  return(enrichment)
}

#' Create SingleCellExperiment object
#'
#' @param expr_data gene x sample expression matrix.
#' @param cell_data sample x k matrix, where k is the number of
#'   cell-dependent features such as treatment, batch number, etc.
#' @param gene_data gene x m matrix, where m is the number of gene features,
#'   such as gene length, SPU number, etc.
#' @param controls optional character vector where elements in the vector
#'   denote control cells by looking at the `treatment` value in the `cell_data`
#'   matrix.
#'
#' @return SingleCellExperiment
#' @export
#'
#' @examples
create_single_cell_object <- function(expr_data, cell_data, gene_data,
                                      controls=c('ASW', 'DMSO')) {
  row.names(gene_data) <- gsub('model', 'TU', row.names(gene_data))
  if (!is.null(controls)) {
    cell_data$controls <- sapply(cell_data$treatment, function(x) {
      as.character(x) %in% controls
      })
  }
  bad_cells <- c("ASW_C07_2016-07-29", "ASW_C01_2016-07-29", "ASW_C09_2016-07-29",
                 "ASW_A02_2016-07-29", "ASW_G06_2016-07-29", "ASW_E05_2016-07-29",
                 "ASW_A12_2016-07-29", "ASW_G03_2016-07-29", "ASW_E04_2016-07-29",
                 "ASW_E11_2016-07-29", "ASW_C02_2016-07-29",
                 "Chlorate-PMCs-1-G05_2018-07-01",
                 "Chlorate-PMCs-1-E06_2018-07-01", "MK886-PMCs-2-A05_2018-07-01",
                 "MK886-PMCs-2-F10_2018-07-01", "MK886-PMCs-3-E05_2018-07-01",
                 "MK886-PMCs-3-D01_2018-07-01", "MK886-PMCs-3-B11_2018-07-01",
                 "MK886-PMCs-3-G10_2018-07-01", "MK886-PMCs-3-G05_2018-07-01")

  keep_genes <- intersect(row.names(expr_data), row.names(gene_data))
  expr_data <- expr_data[keep_genes, ]
  gene_data <- gene_data[keep_genes, ]
  expr_data <- expr_data[!(colnames(expr_data) %in% bad_cells)]
  cell_data <- cell_data[colnames(expr_data), ]
  gene_data$UniProt.Name <- sapply(gene_data$UniProt.Name, format_uniprot_name)
  sce <- SingleCellExperiment::SingleCellExperiment(
                                 assays=list('logcounts'=as.matrix(expr_data)),
                                 colData=as.matrix(cell_data),
                                 rowData=as.matrix(gene_data))
  reducedDim(sce, 'umap') <- umap::umap(t(logcounts(sce)))$layout
  return(sce)
}

format_uniprot_name <- function(uniprot_name) {
  return(strsplit(as.character(uniprot_name), '_')[[1]][1])
}

variable_genes <- function(sce) {
  fit <- trendVar(sce, use.spikes=FALSE)
  decomp <- decomposeVar(sce, fit)
  top.hvgs <- order(decomp$bio, decreasing=TRUE)
  plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
  o <- order(decomp$mean)
  lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
  points(fit$mean, fit$var, col="red", pch=16)
}


#' Plot Gene Expression
#'
#' Plots cells in a reduced dimension, where cells are colored accor$ding to
#' their expression level of a specified gene.
#'
#' @param sce SingleCellExperiment containing expression data, cell data, and
#'   gene data.
#' @param gene_id character argument specifying identifier for gene of
#'   interest (i.e. Notch).
#' @param id_column character argument specifying which gene data column the
#'   identifier is from.
#' @param title Optional character argument to title graph.
#'
#' @return ggplot object with cells plotted along UMAP dimensions.
#' @export
#'
#' @examples
plot_gene_expression <- function(sce, gene_id, id_column, title=NULL) {
  if (is.null(reducedDim(sce, 'umap'))) {
    reducedDim(sce, 'umap') <- umap(t(logcounts(sce)))$layout
  }
  cell_plot <- as.data.frame(reducedDim(sce, 'umap'))
  colnames(cell_plot) <- c('UMAP1', 'UMAP2')
  evm <- row.names(rowData(sce))[which(rowData(sce)[ , id_column] == gene_id)]
  if (length(evm) > 1) {
    warning(paste0('Multiple entries for ', gene_id,
                   '. Using most-highly expressed.'))
    evm_sum <- apply(exprs(sce)[evm, ], 1, sum)
    evm <- evm[order(-evm_sum)]
  }
  cell_plot$Expression <- as.vector(logcounts(sce)[evm, ])
  p <- ggplot(cell_plot, aes(x = UMAP1, y = UMAP2, col = Expression)) +
    geom_point() +
    scale_color_viridis_c(name='log2(tpm)')
  if (! is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

subset_on_gene_list <- function(sce, gene_list, id_column, k=NULL, max_k=15) {
  row_nums <- which(rowData(sce)[ , id_column] %in% gene_list)
  expr_data <- as.data.frame(t(exprs(sce)[row_nums, ]))
  if (is.null(k)) {
    print("Esimating number of clusters...")
    gap_statistic <- as.data.frame(cluster::clusGap(expr_data, cluster::pam,
                                                    K.max=max_k, B=500)$Tab)
    k_estimate <- sapply(1:(nrow(gap_statistic) - 1), function(i) {
      gap_statistic$gap[i] - gap_statistic$gap[i + 1] + gap_statistic$SE.sim[i]
    })
    k <- max_k
    sign_change <- which(k_estimate > 0)
    if (length(sign_change) > 0) {
      k <- sign_change[1]
    }
    # smoothed <- stats::loess.smooth(seq(1:max_k), gap_statistic$gap)
    # idx <- inflection::findiplist(smoothed$x, smoothed$y, 0)[2,2]
    # k <- as.integer(smoothed$x[idx])
  }
  medoids <- cluster::pam(expr_data, k=k)
  subset_umap <- umap::umap(expr_data)$layout
  cell_data <- as.data.frame(colData(sce))
  cell_data$cluster <- as.factor(medoids$clustering)
  cell_data$umap1 <- subset_umap[ , 1]
  cell_data$umap2 <- subset_umap[ , 2]

  out_list <- list()

  umap_plot <- ggplot(cell_data, aes(x=umap1, y=umap2, col=cluster)) +
    geom_point() +
    ggtitle('Cell Clusters')
  out_list$plot <- umap_plot
  if (exists('gap_statistic')) {
    gap_statistic$n.cluster <- 1:max_k
    k_plot <- ggplot(gap_statistic, aes(x=n.cluster, y=gap)) +
      geom_errorbar(aes(ymin=gap - SE.sim, ymax=gap + SE.sim), width=0.1) +
      geom_point() +
      geom_smooth(se=FALSE, weight=0.2) +
      geom_vline(aes(xintercept=k), color='red', linetype='dashed') +
      ggtitle(paste0('Cluster estimation: inflection point at k=', k))
    out_list$plot <- plot_grid(out_list$plot, k_plot)
  }
  out_list$sce <- SingleCellExperiment::SingleCellExperiment(
                                        assays=list('logcounts'=exprs(sce)),
                                        colData=cell_data,
                                        rowData=rowData(sce))
  out_list$model <- medoids
  return(out_list)
}