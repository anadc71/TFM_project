#' myHeatMap
#'
#' Creates a heatmap using a DESeqDataSet class object as input.
#' By default a hierarquic cluster of the mean expression values of the miRNas is done.
#' Optional, a cluster of the most variable miRNAs against samples.
#' @param dds A DESeqDataSet subclass object
#' @param clustering By default clustering by mean expression (set to TRUE). Clustering by variance when set to FALSE.
#' @param colors Define palette for the heatmap.
#' @param ann_colors Define palette for annotations
#' @param mapTitle The title for the heatmap
#' @importFrom DESeq2 rlog
#' @importFrom pheatmap pheatmap
#' @importFrom genefilter rowVars
#' @return None
#' @examples
#' cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#' conditions <- factor(rep(1:2, each=5))
#' dds <- DESeqDataSetFromMatrix(cnts, DataFrame(conditions), ~ conditions)
#' dds<-DESeq(dds)
#' df<- as.data.frame(DataFrame(dds)[,c("conditions","replicate")])
#' df<-data.frame(conditions)
#' colors<-colorRampPalette(c("blue", "white", "green3"))(256)
#' ann_colors<-list(V1=c("#1B9E77","#D95F02"))
#' myHeatMap(dds, clustering="FALSE", colors, ann_colors, mapTitle="Example")
#' @export
myHeatMap<-function(dds,clustering=TRUE, colors, ann_colors, mapTitle){
    rld_dds<-rlog(dds, blind=TRUE)
    if (clustering){
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:50]
    map_data<-assay(rld_dds)[select,]
    pheatmap(map_data, cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df,
             annotation_colors = ann_colors, col=colors,
             fontsize_row=6, fontsize_col=9, main=mapTitle)
  } else {
    topVar<-head(order(rowVars(assay(rld_dds)), decreasing=TRUE),50)
    mat<-assay(rld_dds)[topVar,]
    mat<-mat-rowMeans(mat)
    pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df,
             annotation_colors = ann_colors, col=colors,
             fontsize_row=6, fontsize_col=9, main=mapTitle)
  }
  return(heatmap)
}
