#' myPlotCounts
#'
#' Creates a boxplot of the most significative upregulated or downregulated feature
#'
#' @param dds A DESeqDataSet subclass object
#' @param results A result table from a DESeq analysis
#' @importFrom DESeq2 plotCounts
#' @importFrom ggplot2 ggplot
#' @return None
#' @examples
#' dds <- DESeqDataSetFromMatrix(m=4)
#' dds<-DESeq(dds)
#' res<-results(dds, contrast=c("condition","B","A"))
#' @export
myPlotCounts<-function(dds, results, FoldChange=TRUE){
    if(FoldChange){
      upSig<-subset(res, padj < 0.01 & log2FoldChange > 1.5)
      ordupSig<-upSig[order(upSig$padj),]
    topmiRNA<-rownames(ordupSig)[which.min(ordupSig$padj)]
    plot<-plotCounts(dds, gene=topmiRNA, intgroup="conditions", returnData = TRUE)
    box_plot<-ggplot(plot, aes(x=conditions,y=count, fill=conditions))+geom_boxplot()
    box_plot+ scale_fill_manual(values=c(control="wheat2", treatment="gray51"))+ggtitle (paste0(topmiRNA))+
      theme_classic()+theme(legend.position='none')+
      xlab("Conditions") + ylab("Normalized counts")
  }else{
    dnSig<-subset(res, padj < 0.01 & log2FoldChange < -1.5)
    orddnSig<-dnSig[order(dnSig$padj),]
    topmiRNA<-rownames(orddnSig)[which.min(orddnSig$padj)]
    plot<-plotCounts(dds, gene=topmiRNA, intgroup="conditions", returnData = TRUE)
    box_plot<-ggplot(plot, aes(x=conditions,y=count, fill=conditions))+geom_boxplot()
    box_plot+ scale_fill_manual(values=c(control="wheat2", treatment="gray51"))+ggtitle (paste0(topmiRNA))+
      theme_classic()+theme(legend.position='none')+
      xlab("Conditions") + ylab("Normalized counts")
  }
}


