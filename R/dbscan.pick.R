#' pick internal spike-in genes
#'
#' @description pick internal spike-in genes, using \code{\link{cluster.count}}.
#' @param dis gene distance from \code{\link{calculate.dis}}.
#' @param ngene a series of expected number of IS genes.
#' @import dbscan
#' @param solution specify the increasing rate of scanning radius. Default 100.
#' @export
#' @return internal spike-in genes picked up
#' @seealso \code{\link{cluster.count}} and \code{\link{calculate.dis}}
#' @examples
#' ## Not run
#' spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5,solution=100)
#' ## End(Not run)

dbscan.pick<-function(dis,ngene=(1:floor(nrow(dis)/25))*5,solution=100){
  genelist<-rownames(dis)
  ngene<-ngene[ngene<=length(genelist)]
  ngene<-sort(ngene)
  dis<-as.dist(dis)
  xx<-as.vector(dis)
  xx<-xx[xx!=0]
  step<-min(xx)
  ngene_index<-1
  output<-list()
  while(ngene_index<=length(ngene)){
    cluster<-0
    while(max(cluster.count(cluster))<ngene[ngene_index]){
      step<-step*10
      cluster<-dbscan(dis,eps=step)$cluster
    }
    step<-step/10
    cluster<-0
    i<-1
    while(i<=solution&ngene_index<=length(ngene)){
      i<-i+1
      cluster<-dbscan(dis,eps=step+step*9/solution*i)$cluster
      if(max(cluster.count(cluster))<ngene[ngene_index]) next
      ngene_index2<-max(which(max(cluster.count(cluster))>=ngene))
      while(ngene_index<=ngene_index2){
        output[[ngene_index]]<-genelist[cluster==which.max(cluster.count(cluster))]
        ngene_index<-ngene_index+1
      }
    }
  }
  output<-output[c(T,!sapply(2:length(output),function(x) identical(output[[x-1]],output[[x]])))]
  return(output)
}  ##pick internal spike-in genes
