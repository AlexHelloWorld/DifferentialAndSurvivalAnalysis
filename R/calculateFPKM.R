#' calculate FPKM from countTable and coldata
#' @return FPKM matrix
#' @param countTable a raw read count table
#' @param coldata an annotation table for samples, must include attributes naming "barcode", "patient_barcode", "survival", "vital_status", "total_count", "cohort"
#' @param geneInfo an annotation table for genes, must include attribute nameing "gene_length"
#' @export
calculateFPKMMatrix <- function(countTable, coldata, geneInfo){
  for(i in 1:ncol(countTable)){
    totalCount <- coldata[colnames(countTable)[i], "total_count"]
    countTable[, i] <- countTable[, i]*10^9/totalCount
  }
  for(i in 1:nrow(countTable)){
    geneLength <- geneInfo[geneInfo$gene == rownames(countTable)[i], "gene_length"]
    countTable[i, ] <- countTable[i, ]/geneLength
  }
  countTable
}

#' plot RPKM boxplot
#' @return boxplot of FPKM by ggplot2
#' @param FPKMTable FPKM matrix
#' @param coldata an annotation table for samples, must include attributes naming "barcode", "patient_barcode", "survival", "vital_status", "total_count", "cohort"
#' @param geneInfo an annotation table for genes, must include attribute nameing "gene_length"
#' @param attribute attribute used to group data (x axis)
#' @param geneName the name of gene to draw boxplot
#' @param geneId the Id of gene to draw boxplot
#' @import ggplot2
#' @import ggpubr
#' @export
plotFPKMBoxplot <- function(FPKMTable, coldata, geneInfo, attribute, geneName = "", geneId = ""){
  if(geneId == "" & geneName == ""){
    print("Must provide either of geneId or geneName")
    return
  }else if(geneId != ""){
    geneName <- as.character(geneInfo[geneInfo$gene == geneId, "gene_name"])
  }else{
    geneId <- geneInfo[geneInfo$gene_name == geneName, "gene"]
  }
  fpkm <- data.frame(FPKM = unlist(FPKMTable[as.character(geneId), ]), id = colnames(FPKMTable))
  coldata$id <- rownames(coldata)
  fpkm <- merge(fpkm, coldata[, c("id", attribute)], by = "id")
  fpkm <- fpkm[!is.na(fpkm[, attribute]), ]

  p <- ggplot(fpkm, aes(x = as.factor(fpkm[[attribute]]), y=FPKM, color = as.factor(fpkm[[attribute]]))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.1, height = 0)) +
    ggtitle(paste(geneId, geneName, sep = "\t")) +
    labs(x = attribute, color = attribute)

  p + stat_compare_means()
}
