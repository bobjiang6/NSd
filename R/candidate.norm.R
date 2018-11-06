#' normalization by spike candidate
#'
#' @description calculate pairwise DoR distance between genes, using \code{\link{calculate.dor}}
#' @param mat the expression matrix.
#' @param spike_candidate output from \code{\link{dbscan.pick}}
#' @param ncore the number of cores used. Default 1.
#' @import parallel
#' @export
#' @return  a list containing the normalization results for each candidate set
#' @seealso candidate_res
#' @examples
#' ## Not run
#' candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=4)
#' ## End(Not run)

candidate.norm<-function(mat,spike_candidate,ncore=1){
  cl <- makeCluster(ncore)
  candidate_ref<-parLapply(cl,spike_candidate,estimate.ref,mat=mat)
  temp<-lapply(candidate_ref,function(x) apply(mat,2,estimate.sf,ref=x))
  stopCluster(cl)
  candidate_res<-list(sf=sapply(temp,function(x) x[1,]),inst=sqrt(sapply(temp,function(x) x[2,])),ngene=sapply(temp,function(x) x[3,]))
  candidate_res$spike<-spike_candidate

  return(candidate_res)
}  ##normalization by spike candidate
