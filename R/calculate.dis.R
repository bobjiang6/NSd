#' calculate pairwise DoR distance between genes
#'
#' @description calculate pairwise DoR distance between genes, using \code{\link{calculate.dor}}
#' @param mat the expression matrix.
#' @param dectection_rate the rate set for dectecing distance. If the running time is too long, you can try detection_rate=0.95 or other higher values to filter more genes. Default 0.9.
#' @param ncore the number of cores used. Default 1.
#' @export
#' @return an matrix with pairwise DoR distance
#' @seealso \code{\link{calculate.dor}}
#' @examples
#' ## Not run
#' mat<-as.matrix(read.csv(file="GSM1599494_ES_d0_main.csv",sep=",",header=F,row.names=1))
#' gene_dis<-calculate.dis(mat=mat,detection_rate=0.9,ncore=4)
#' ## End(Not run)

calculate.dis<-function(mat,detection_rate=0.9,ncore=1){
  mat2<-mat[apply(mat,1,function(x) sum(x>0)/length(x))>detection_rate,]
  gene_dis<-matrix(0,nrow(mat2),nrow(mat2))
  cl <- makePSOCKcluster(ncore)

  for(i in 2:nrow(mat2)-1){
    gene_dis[i,i:nrow(mat2)]<-parallel::parRapply(cl=cl,mat2[i:nrow(mat2),],calculate.dor,yy=mat2[i,])
    gene_dis[i:nrow(mat2),i]<-gene_dis[i,i:nrow(mat2)]
  }
  stopCluster(cl)
  rownames(gene_dis)<-rownames(mat2)
  colnames(gene_dis)<-rownames(mat2)
  return(gene_dis)
}
