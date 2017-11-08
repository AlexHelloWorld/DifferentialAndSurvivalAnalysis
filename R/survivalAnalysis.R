#' Prepare table for survival analysis
prepareTableForSurvivalAnalysis <- function(gene, rpkm, coldata){
  #clean rpkm table
  expression <- data.frame(t(rpkm[gene, ]))
  colnames(expression) <- "rpkm"
  expression$barcode <- rownames(expression)

  #clean up coldata
  clinical <- subset(coldata, select=c("vital_status", "survival"))
  clinical$barcode <- rownames(clinical)
  clinical <- subset(clinical, subset = (!is.na(clinical$survival) & clinical$survival>0))
  clinical$survivalCensor <- NA
  for(i in 1:nrow(clinical)){
    if(clinical$vital_status[i] == "Dead"){
      clinical$survivalCensor[i] <- 1
    }else{
      clinical$survivalCensor[i] <- 0
    }
  }
  clinical <- merge(clinical, expression)
  clinical
}

#' Cox survival analysis
#' @return a list containing cox regression p value and coef
#' @param gene gene id
#' @param rpkm FPKM matrix
#' @param coldata an annotation table for samples
#' @import survival
#' @export
coxSurvivalAnalysis <- function(gene, rpkm, coldata){
  clinical <- prepareTableForSurvivalAnalysis(gene, rpkm, coldata)
  #do cox survival analysis
  survFit <- coxph(Surv(survival, survivalCensor)~rpkm, data = clinical)
  survSummary <- summary(survFit)
  if(length(survSummary$waldtest) != 3){
    list(pValueCox = NA, coefCox = NA)
  }else{
    list(pValueCox = survSummary$waldtest[[3]], coefCox = survSummary$coefficients[[1]])
  }
}

#' Cox survival analysis for a list of genes
#' @return a data frame containing cox regression results (gene id, cox regression p value, coef)
#' @param geneList a list of gene id
#' @param rpkm FPKM matrix
#' @param coldata coldata an annotation table for samples
#' @export
coxSurvivalAnalysisForGenes <- function(geneList, rpkm, coldata){
  firstResult <- coxSurvivalAnalysis(geneList[1], rpkm, coldata)
  output <- data.frame(gene = geneList[1], pValueCox = firstResult$pValueCox, coefCox = firstResult$coefCox)
  if(length(geneList) >=2 ){
    for(i in 2:length(geneList)){
      result <- coxSurvivalAnalysis(geneList[i], rpkm, coldata)
      output <- rbind(output, data.frame(gene = geneList[i], pValueCox = result$pValueCox, coefCox = result$coefCox))
    }
  }
  output
}

#' plot survival curves
#' @return a list containing survival plot, p value and table for plotting
#' @param gene id
#' @param coldata an annotation table for samples
#' @param topn group samples into two groups: high expression and low expression. Classify the top n samples (by FPKM) into high expression group.
#' @param rpkmCutoff classify samples into two groups by a rpkmCutoff
#' @import survminer
#' @import survival
#' @export
plotSurvivalAnalysis <- function(gene, rpkm, coldata, topn = 10, rpkmCutoff = 0){
  clinical <- prepareTableForSurvivalAnalysis(gene, rpkm, coldata)
  if(rpkmCutoff == 0){
    clinical <- clinical[order(clinical$rpkm, decreasing = TRUE), ]
    clinical$high_expression <- c(rep(TRUE, topn), rep(FALSE, nrow(clinical)-topn))
  }else{
    clinical <- clinical[order(clinical$rpkm, decreasing = TRUE), ]
    clinical$high_expression <- (clinical$rpkm > rpkmCutoff)
  }
  survFit <- survfit(Surv(survival, survivalCensor)~high_expression, data=clinical)
  plot <- ggsurvplot(survFit, data=clinical, pval = TRUE)
  survSummary <- summary(survFit)
  survDiff <- survdiff(Surv(survival, survivalCensor)~high_expression, data=clinical)
  list(survPlot = plot, survDiff = survDiff, plotTable = clinical)
}

