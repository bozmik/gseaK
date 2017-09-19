#' Gene Set
#'
#' GSet
#'
#' @param gSets gene Set
#' @param gene_name_ord ordered gene names
#' @param ord ordered values of statistic
#' @param n_perm number of permutations
#' @param expr expression dataset
#' @param class phenotype data
#' @param index_ph shuffled phenotype labels
#' @return returns a \code{vector} object
#' @useDynLib gseaK
#' @importFrom Rcpp sourceCpp
#' @docType methods
#' @rdname GSet-methods
#' @export
setGeneric("GSet",
           function(gSets, gene_name_ord, ord, n_perm, expr, class, index_ph)
             standardGeneric("GSet") )

#' @aliases GSet-method
#' @rdname GSet-methods
setMethod("GSet", signature("vector","character", "vector","numeric","matrix","factor","matrix"),
          function(gSets, gene_name_ord, ord, n_perm, expr, class, index_ph){

result<-matrix(nrow = 1,ncol=8)
da=gSets[!is.na(gSets)]  ##Gene set
#name_gene<-colnames(gSets)  ##name of GS
#splitt<-strsplit(name_gene,'path.')
#name_gene_set<-splitt[[1]][2]
#colnames(da)<-name_gene_set

N_H<-length(da)  #n in GS
N<-length(gene_name_ord)  #N ordered genes

#library(st)
#library(Rcpp)
#sourceCpp("ginGS.cpp")
#####Genes in GS
l<-ginGS(gene_name_ord, as.character(da), ord);

N_R<-sum(abs(as.numeric((l$stat))))
P_miss<-(1/(N-N_H))

#sourceCpp("ES.cpp")
ES_matrix<- ES(ord, gene_name_ord, P_miss, N_R, l$poz)
colnames(ES_matrix)=c("P_hit","P_miss","ES")
rownames(ES_matrix)<-gene_name_ord

sum_ES=max(abs(ES_matrix[,3]),na.rm = T)
x_ES<-which(abs(ES_matrix[,3])==sum_ES)

ES_obs<-ES_matrix[x_ES,3]

######...plots

ES_p=matrix(nrow=n_perm)

###permutations

#library(doParallel)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
ES_p <- foreach(i=1:n_perm, .combine=c, .noexport = c("ginGS","ES")) %dopar% {
#  library(st)
 # library(Rcpp)
  expr2 <- expr[,index_ph[i,]]
  p_stat_t_p <- as.matrix(diffmean.stat(t(expr2),class))

  p_stat_t_p <- as.data.frame(abs(p_stat_t_p))

  ord_p <- p_stat_t_p[order(p_stat_t_p,decreasing = T),,drop=F]   ###ranking of genes
  gene_name_ord_p<-as.matrix(rownames(ord_p))

  N<-dim(ord_p)[1]

  #####Genes in GS

 # sourceCpp("ginGS.cpp")

  l_p<-ginGS(gene_name_ord_p, as.character(da), ord_p[,1]);

  N_R_p<-sum(abs(as.numeric((l_p$stat))))

#  sourceCpp("ES.cpp")

  ES_matrix_p<-ES(ord_p[,1], rownames(ord_p), P_miss,N_R_p,l_p$poz)

  rownames(ES_matrix_p)<-gene_name_ord_p

  max_ES_p <- max(abs(ES_matrix_p[,3]),na.rm = T)
  x_ES_p<-which(abs(ES_matrix_p[,3])==max_ES_p)
  ES_matrix_p[x_ES_p,3]  ## obserwowane ES z kazdej permtacji
}
stopCluster(cl)


###plots


### p-value
#library(sfsmisc)
dens <- density(ES_p)
##### poprawić!!!
if(ES_obs>=0){
  p<-length(which(ES_p>=ES_obs))
  p_val<-p/n_perm
  if (mean(abs(ES_p[ES_p>=0]))=='NaN') NES=0 else NES<-ES_obs/mean(ES_p[ES_p>=0])
  lll<-which(dens$x>ES_obs)
  if (length(lll)==0) p_val_kern<-0 else p_val_kern<-integrate.xy(dens$x[lll],dens$y[lll])
}else{
  p<-length(which(ES_p<=ES_obs))
  p_val<-p/n_perm
  if (mean(abs(ES_p[ES_p<0]))=='NaN') NES=0 else NES<-ES_obs/mean(abs(ES_p[ES_p<0]))
  lll<-which(dens$x<ES_obs)
  if (length(lll)==0) p_val_kern<-0 else p_val_kern<-integrate.xy(dens$x[lll],dens$y[lll])
}


## NES_p<-fun(ES_p)
##  NES_p<-cbind(ES_p,NES_p)
##    colnames(NES_p)<-c('ES_p','NES_p')


##   if(NES>=0){
#    fr<-length(which(NES_p[,2]>=NES))  #fraction
#    fr_poz<-length(which(NES_p[,2]>=0))
#    fraction=fr/fr_poz
#  }else{
#    fr<-length(which(NES_p[,2]<=NES))
#    fr_neg<-length(which(NES_p[,2]<0))
#    fraction=fr/fr_neg
#  }

#result[1]<-name_gene_set
result[2]<-ES_obs
result[3]<-p_val
#    result[t,4]<-NES ##NES observed
##  result[t,6]<-fraction
result[7]<-p_val_kern

return(result)

          }
)
