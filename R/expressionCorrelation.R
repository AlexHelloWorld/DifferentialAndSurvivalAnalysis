#' calculate correlation of expression between two genes
expressionCorrelationSingleGene <- function(gene1, matrix1, gene2, matrix2){
  #presume that matrix1 and matrix2 have the same rownames
  gene1Expression <- unlist(matrix1[gene1, ])
  gene2Expression <- unlist(matrix2[gene2, ])
  
  corTest <- cor.test(gene1Expression, gene2Expression, method = "pearson")
  c(corTest$estimate, corTest$p.value)
}

#' calculate correlation of expression between a gene and a gene list
#' @return a data.frame containing gene id, pearson coefficient, p value
#' @param gene targe gene id
#' @param matrix1 FPKM matrix containing target gene id
#' @param matrix2 FPKM matrix for looping and calculate pearson correlation with target gene
#' @export
expressionCorrelation <- function(gene, matrix1, matrix2){
  #make sure two matrices have the same columns
  names <- intersect(colnames(matrix1), colnames(matrix2))
  matrix1 <- matrix1[, names]
  matrix2 <- matrix2[, names]
  firstResult <- expressionCorrelationSingleGene(gene, matrix1, rownames(matrix2)[1], matrix2)
  output <- data.frame(gene = rownames(matrix2)[1], cor = firstResult[1], pValue = firstResult[2])
  if(nrow(matrix2) >= 2){
    for(i in 2:nrow(matrix2)){
      result <- expressionCorrelationSingleGene(gene, matrix1[gene, ], rownames(matrix2)[i], matrix2[i, ])
      output <- rbind(output, data.frame(gene = rownames(matrix2)[i], cor = result[1], pValue = result[2]))
    }
  }
  output
}
