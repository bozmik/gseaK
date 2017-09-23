#' Gene Set
#'
#' GSet
#'
#' @param gSets gene Set
#' @param gene_name_ord ordered gene names
#' @param ord ordered values of statistic
#' @param n.perm number of permutations
#' @param expr expression dataset
#' @param class phenotype data
#' @param index_ph shuffled phenotype labels
#' @param stest Rank metrics
#' @param abs absolute value
#' @param plot1 plot1
#' @param plot2 plot2
#' @param plot3 plot3
#' @return returns a \code{vector} object
#' @useDynLib gseaK
#' @importFrom Rcpp sourceCpp
#' @docType methods
#' @rdname GSet-methods
#' @export
setGeneric("GSet",
           function(gSets, gene_name_ord, ord, n.perm, expr, class, index_ph, stest, abs, plot1, plot2, plot3)
             standardGeneric("GSet") )

#' @aliases GSet-method
#' @rdname GSet-methods
setMethod("GSet", signature("vector","character", "vector"),
          function(gSets, gene_name_ord, ord, n.perm, expr, class, index_ph, stest, abs, plot1, plot2, plot3){

name_gene <- gSets[1]
gSets <- as.numeric(gSets[-1])
result<-matrix(nrow = 1,ncol=7)
da=as.numeric(gSets[!is.na(gSets)])  ##Gene set
#name_gene<-colnames(gSets)  ##name of GS
splitt<-strsplit(name_gene,'path.')
name_gene_set<-splitt[[1]][2]
#colnames(da)<-name_gene_set
N_H<-length(da)  #n in GS
N<-length(gene_name_ord)  #N ordered genes

#library(st)
#library(Rcpp)
#sourceCpp("ginGS.cpp")
#####Genes in GS
l<-ginGS(gene_name_ord,as.character(da),ord)


N_R<-sum(abs(as.numeric((l$stat))))
P_miss<-(1/(N-N_H))

ES_matrix<- ES(ord, gene_name_ord, P_miss, N_R, l$pos)
colnames(ES_matrix)=c("P_hit","P_miss","ES")
rownames(ES_matrix)<-gene_name_ord

if(plot1==TRUE){
  plotDis(ES_matrix, name_gene_set)}

sum_ES=max(abs(ES_matrix[,3]),na.rm = T)
x_ES<-which(abs(ES_matrix[,3])==sum_ES)

ES_obs<-ES_matrix[x_ES,3]

if(plot2 == TRUE){
  plotES(ES_matrix, ES_obs, x_ES, l$pos, name_gene_set)}

ES_p=matrix(nrow=n.perm)

###permutations

cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
ES_p <- foreach(i=1:n.perm, .combine=c, .noexport = c("ginGS","ES")) %dopar% {

  expr2 <- expr[,index_ph[i,]]

  ###### Rank metrics
  if(stest == "FC"){
    p_stat_t_p <- as.matrix(diffmean.stat(t(expr2),class))
  } else if(stest == "shrinkage.t"){
    p_stat_t_p <- as.matrix(shrinkt.stat(t(expr2),class,var.equal=F,paired=FALSE, verbose=T))
  } else if(stest == "2s.Bayesian"){
    p_stat_t_p <- as.matrix(RankingFoxDimmic(expr2,class,type = "unpaired")@statistic)
  } else if(stest == "Efron"){
    p_stat_t_p <- as.matrix(RankingEbam(expr2,class,type = "unpaired")@statistic)
  } else if(stest == "SAM"){
    p_stat_t_p <- as.matrix(RankingSam(expr2,class,type = "unpaired")@statistic)
  } else if(stest == "penalized.t"){
    p_stat_t_p <- as.matrix(RankingSoftthresholdT(expr2,class, type = c("unpaired"))@statistic)
  } else if(stest == "bayesian.t"){
    p_stat_t_p <- as.matrix(RankingBaldiLong(expr2,class, type = c("unpaired"))@statistic)
  } else if(stest == "moderatet.t"){
    p_stat_t_p <- as.matrix(RankingLimma(expr2,class,type = c("unpaired"))@statistic)
  }
  ##else if(stest == "moderated.wt"){
  #  p_stat_t_p <- as.matrix(mwt(expr2,class,log.it = FALSE)$MWT)
  #}
  rownames(p_stat_t_p)=rownames(expr2)

  if (abs==TRUE) {
    p_stat_t_p<-as.data.frame(abs(p_stat_t_p))
  } else {
    p_stat_t_p<-as.data.frame(p_stat_t_p)}

  ord_p <- p_stat_t_p[order(p_stat_t_p,decreasing = T),,drop=F]   ###ranking of genes
  gene_name_ord_p<-as.matrix(rownames(ord_p))

  N<-dim(ord_p)[1]

  #####Genes in GS

  l_p<-ginGS(gene_name_ord_p, as.character(da), ord_p[,1]);

  N_R_p<-sum(abs(as.numeric((l_p$stat))))

  ES_matrix_p<-ES(ord_p[,1], rownames(ord_p), P_miss,N_R_p,l_p$pos)

  rownames(ES_matrix_p)<-gene_name_ord_p

  max_ES_p <- max(abs(ES_matrix_p[,3]),na.rm = T)
  x_ES_p<-which(abs(ES_matrix_p[,3])==max_ES_p)
  ES_matrix_p[x_ES_p,3]  ## obserwowane ES z kazdej permtacji
}
stopCluster(cl)

dens <- density(ES_p)

if(plot3 == TRUE){
  plotPerm(ES_p, ES_obs, dens, name_gene_set)}


### p-value

if(ES_obs>=0){
  p<-length(which(ES_p>=ES_obs))
  p_val<-p/n.perm
  if (mean(abs(ES_p[ES_p>=0]))=='NaN') NES=0 else NES<-ES_obs/mean(ES_p[ES_p>=0])
  lll<-which(dens$x>ES_obs)
  if (length(lll)==0) p_val_kern<-0 else p_val_kern<-integrate.xy(dens$x[lll],dens$y[lll])
}else{
  p<-length(which(ES_p<=ES_obs))
  p_val<-p/n.perm
  if (mean(abs(ES_p[ES_p<0]))=='NaN') NES=0 else NES<-ES_obs/mean(abs(ES_p[ES_p<0]))
  lll<-which(dens$x<ES_obs)
  if (length(lll)==0) p_val_kern<-0 else p_val_kern<-integrate.xy(dens$x[lll],dens$y[lll])
}

mp<-mean(ES_p[ES_p>=0])
mn<-mean(abs(ES_p[ES_p<0]))

NES_p <- NES(ES_p, mp, mn)

if(NES>=0){
  fr<-length(which(NES_p >= NES))  #fraction
  fr_poz<-length(which(NES_p >= 0))
  fraction=fr/fr_poz
}else{
  fr<-length(which(NES_p <= NES))
  fr_neg<-length(which(NES_p < 0))
  fraction=fr/fr_neg
}

result[1]<-ES_obs
result[2]<-p_val
result[3]<-NES ##NES observed
result[5]<-fraction
result[6]<-p_val_kern

return(result)

  }
)
