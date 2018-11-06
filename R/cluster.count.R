#' count the size of a cluster
#'
#' @description count cluster size, used for \code{\link{dbscan.pick}}
#' @param x a cluster for count
#' @export
#' @return  size of the cluster
#' @seealso \code{\link{dbscan.pick}}
#' @examples
#' ## Not run
#' ##TBD
#' ## End(Not run)
cluster.count<-function(x){
  output<-numeric()
  for(i in 1:max(x)){
    output[i]<-sum(x==i)
  }
  return(output)
}  ##count cluster size
