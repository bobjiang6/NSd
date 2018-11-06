#' calculate DoR between two vector
#'
#' @description calculate DoR between two vector
#' @param xx
#' @param yy
#' @export
calculate.dor<-function(xx,yy){
  retain<-xx!=0&yy!=0
  ratio<-xx[retain]/yy[retain]
  output<-sqrt(sum(log2(ratio/median(ratio))^2)/(sum(retain)-1))
  return(output)
}
