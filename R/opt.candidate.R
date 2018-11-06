#' Choose an appropriate IS geneset
#'
#' @description Choose an appropriate IS geneset based on the results from candidate_res.
#' @param mat the expression matrix.
#' @param candidate_res output from \code{\link{candidate.norm}}
#' @param baseline_threshold specifies the threshold of instability score used to choose the baseline geneset (see our article for more details).
#' @param p_value specifies the threshold of p value for F-test. Default 0.05.
#' @export
#' @return  a list with 5 elements.The output ISnorm_res$normalized contains the normalized matrix. The output ISnorm_res$size_factor contains the size factor for each cell. These two outputs can be used as inputs for other scRNA-seq analysis. The output ISnorm_res$ISgenes contains the name of IS genes used for normalization. The output ISnorm_res$inst_cell contains the instability score for each cell. The output ISnorm_res$picked contains the index of the optimized candidate geneset. A message will also tell you the index of the chosen geneset.
#' @seealso \code{\link{candidate.norm}}
#' @examples
#' ## Not run
#' candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=4)
#' ## End(Not run)

opt.candidate<-function(mat,candidate_res,baseline_threshold=0.1,p_value=0.05){
  instability<-apply(candidate_res$inst,2,mean)
  if(sum(instability<baseline_threshold)==0){
    baseset<-1
  }
  else{
    baseset<-max((1:length(candidate_res$spike))[instability<0.1])
  }
  inst<-candidate_res$inst
  ngene<-candidate_res$ngene
  pmat<-sapply(baseset:ncol(inst),function(x) pf((inst[,baseset]/inst[,x])^2,df1=ngene[,baseset]-1,df2=ngene[,x]-1))
  picked<-max(which(apply(pmat,2,function(x) (sum(x<p_value)/length(x))<p_value)))+baseset-1
  cat("Candidate set",picked,"is chosen.\n")

  expr<-sweep(mat,2,candidate_res$sf[,picked],FUN="/")
  inst_cell<-candidate_res$inst[,picked]
  ISgenes<-candidate_res$spike[[picked]]
  return(list(normalized=expr,size_factor=candidate_res$sf[,picked],ISgenes=ISgenes,inst_cell=inst_cell,picked=picked))
}
