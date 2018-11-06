#' calculate reference expression given internal spike-in genes
#'
#' @description calculate reference expression given internal spike-in genes
#' @param mat the expression matrix.
#' @param ref_gene the location of spike-in genes in the dataset mat
#' @export
#' @return estimated expression data fot spike-in genes
#' @examples
#' ## Not run
#' ## TBD
#' ## End(Not run)

estimate.ref<-function(mat,ref_gene){
  mat<-mat[ref_gene,]
  cell_sum<-colSums(mat)
  mat<-mat[,cell_sum!=0]
  cell_sum<-cell_sum[cell_sum!=0]
  ref_expr<-apply(sweep(mat,2,cell_sum/median(cell_sum),FUN="/"),1,function(x) mean(x[x!=0]))
  ref_expr<-sort(ref_expr)
  return(ref_expr)
}  ##calculate reference expression given internal spike-in genes
