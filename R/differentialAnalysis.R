#' differential analysis function based on DESeq2
#' Do differential analysis with DESeq2 on single attribute and two attribute
#' @return a \code{\link{DESeqDataSet}} object and output two csv files for differential analysis results and vsd normalized expression table directly to current directory
#' @param countTable a raw read count table
#' @param coldata an annotation table for samples, must include attributes naming "barcode", "patient_barcode", "survival", "vital_status", "total_count", "cohort"
#' @param attribute name of the attribute for differential expression analysis
#' @param geneInfo an annotation table for genes, must include attribute nameing "gene_length"
#' @param datasetIdentifier a prefix for the output csv files
#' @import DESeq2
differentialAnalysisSingleAttribute <- function(countTable, coldata, attribute, geneInfo, datasetIdentifier){
  colnames(coldata) <- replace(colnames(coldata), colnames(coldata) == attribute, "attribute")
  dds <- DESeqDataSetFromMatrix(countData = countTable, colData = coldata, design = ~attribute)
  dds <- DESeq(dds)
  results <- results(dds, alpha=0.05)
  differentialData <- data.frame(gene=row.names(results),
                                 pvalue = results$padj,
                                 lfc = results$log2FoldChange)

  vsdData <- assay(varianceStabilizingTransformation(dds))
  if(length(unique(coldata[["cohort"]])) > 1){
    vsdData <- removeBatchEffect(vsdData, batch = coldata[["cohort"]])
  }
  attributeValues <- unique(coldata[["attribute"]])
  label = paste(datasetIdentifier, attribute, attributeValues[1], attributeValues[2], sep = "_")
  differentialLabel <- paste(label, "differential.csv", sep = "_")
  vsdLabel <- paste(label, "vsd.csv", sep = "_")
  write.csv(merge(differentialData, geneInfo), file = differentialLabel)
  write.csv(vsdData, file = vsdLabel)
  dds
}


#' differential analysis function based on DESeq2
#'
#' Do differential analysis with DESeq2 recursively with attributes defined in 'fields' parameter
#' @return a list of \code{\link{DESeqDataSet}} object and output csv files for differential analysis results and vsd normalized expression table directly to current directory
#' @param countTable a raw read count table
#' @param coldata an annotation table for samples, must include attributes naming "barcode", "patient_barcode", "survival", "vital_status", "total_count", "cohort"
#' @param fields name list of the attributes for differential expression analysis
#' @param geneInfo an annotation table for genes, must include attribute nameing "gene_length"
#' @param datasetIdentifier a prefix for the output csv files
#' @export

differentialAnalysis <- function(countTable, coldata, fields, geneInfo, datasetIdentifier = ""){
  DESeq2List <- list()
  for(attribute in fields){
    inputColdata <- coldata[!is.na(coldata[[attribute]]), ]
    attributeValues <- sort(unique(inputColdata[[attribute]]))
    for(i in 1:length(attributeValues)){
      if(i == length(attributeValues)){
        break()
      }
      for(j in (i+1):length(attributeValues)){
        select <- (inputColdata[[attribute]] == attributeValues[i] | inputColdata[[attribute]] == attributeValues[j])
        selectColdata <- inputColdata[select, ]
        inputCount <- countTable[, rownames(selectColdata)]
        label <- paste(as.character(attribute), as.character(attributeValues[i]),  as.character(attributeValues[j]), sep = "_")
        DESeq2List[[label]] <- differentialAnalysisSingleAttribute(inputCount, selectColdata, attribute, geneInfo, datasetIdentifier)
      }
    }
  }
  DESeq2List
}

#' plot heatmap with vsdData, differentialData
#' @return heatmap by pheatmap
#' @import pheatmap
#' @import limma
#' @import RColorBrewer
#' @export
plotHeatMap <- function(vsdData, coldata, orderAttribute, labelAttribute, differentialData = NULL, pValueCutoff = 1, averageCutoff = 0, varCutoff = 0, clusterRow = TRUE, clusterCol = TRUE, colorBias = 1, showNames = FALSE, rowFontSize = 2, colFontSize = 3){
  #select genes based on significance (pvalue) from differentialData
  if(!is.null(differentialData)){
    rownames(differentialData) <- differentialData$gene
    differentialData <- differentialData[rownames(vsdData), ]
    pSelect <- (!is.na(differentialData$pvalue) & differentialData$pvalue <= pValueCutoff)
    vsdData <- vsdData[pSelect, ]
  }
  #select genes based on FPKM average and variance
  averageSelect <- (rowMeans(vsdData) >= averageCutoff)
  varSelect <- (apply(vsdData, 1, var) >= varCutoff)

  vsdData <- vsdData[averageSelect & varSelect, ]

  #select samples based on attribute
  coldata <- subset(coldata, select = unique(c(orderAttribute, labelAttribute)))
  coldata <- subset(coldata, subset = !is.na(coldata[[orderAttribute]]))
  vsdData <- vsdData[, rownames(coldata)]

  coldata <- coldata[colnames(vsdData), ]
  order <- order(coldata[[orderAttribute]])

  pheatmap(vsdData[, order],
           color = colorRampPalette(rev(brewer.pal(n=10, name="RdYlBu")), bias = colorBias)(100),
           cluster_cols = clusterCol,
           cluster_rows = clusterRow,
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D",
           show_rownames = T,
           show_colnames = showNames,
           fontsize_col = colFontSize,
           fontsize_row = colFontSize,
           annotation_col = coldata)

}
