#' calculate size factor for one cell
#'
#' @description calculate size factor for one cell
#' @param cell cell data
#' @param ref_expr a reference status for cell data
#' @export
#' @return estimated size factor for one cell
#' @examples
#' ## Not run
#' ## TBD
#' ## End(Not run)


estimate.sf<-function(cell,ref_expr){
  cell<-cell[names(ref_expr)]
  sf<-cell[cell!=0]/ref_expr[cell!=0]
  ngene<-length(sf)
  sf_es<-median(sf)
  var<-sum(log10(sf/sf_es)^2)/(length(sf)-1)
  output<-c(sf_es,var,length(sf))
  names(output)<-c("size_factor","var","ngene")
  return(output)
}  ##calculate size factor for one cell
