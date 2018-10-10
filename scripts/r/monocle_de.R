library(monocle)

setwd('PythonModules/scPipe/scripts/python/')

format_uniprot_name <- function(uniprot_name) {
  return(strsplit(as.character(uniprot_name), '_')[[1]][1])
}

main <- function() {
  # read in sample data
  cell_data <- read.csv('louvain_clusters.csv', row.names=1,
                        check.names=FALSE)
  cell_data <- new('AnnotatedDataFrame', data=cell_data)

  # read in expression data -- TPM for monocle
  expr_data <- read.csv('../../output/final/normalized_tpm_matrix.csv',
                        row.names=1, check.names=FALSE)
  filtered_data <- read.csv('variable_genes.csv', row.names=1,
                            check.names=FALSE)
  expr_data <- expr_data[colnames(filtered_data), row.names(filtered_data)]

  # tpm log2(x + 1) transformed, invert
  expr_data <- 2**as.matrix(expr_data) - 1

  # create CellDataSet for monocle analysis
  cell_dataset <- monocle::newCellDataSet(cellData=as.matrix(expr_data),
                                          phenoData=cell_data)
  # get absolute reads for diff test
  rpc_matrix <- monocle::relative2abs(cell_dataset, method='tpm_fraction')

  # need to re-create dataset object with new matrix

  cell_dataset <- monocle::newCellDataSet(cellData=as.matrix(rpc_matrix),
                                          phenoData=cell_data,
                                          lowerDetectionLimit=0.5,
                                          expressionFamily=negbionomial.size())
  diff_test_res <- monocle::differentialGeneTest(cell_dataset,
                                                 fullModelFormulaStr =
                                                   '~louvain',
                                                 relative_expr=FALSE)
  louvain_sig_genes <- subset(diff_test_res, qval < 0.1)
  test <- estimateSizeFactors(cell_dataset)
  test_res <- monocle::differentialGeneTest(test, fullModelFormulaStr = '~louvain')
  test_sig <- subset(test_res, qval < 0.1)
}


