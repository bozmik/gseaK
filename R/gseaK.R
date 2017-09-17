#' Gene Set Enrichment Analysis
#'
#' GSEA
#'
#' @param expr expression data
#' @param pheno phenotype
#' @param gSets Gene Sets
#' @param stest statistical test
#' @param abs absolute value
#' @param kernel Kernal density funcion
#' @param n.perm number of permutations
#' @param correction multiple hypothesis testing
#' @return returns a \code{matrix} object
#' @useDynLib gseaK
#' @importFrom Rcpp sourceCpp
#' @docType methods
#' @rdname gseaK-methods
#' @export
setGeneric("gseaK",
           function(expr, pheno, gSets,
                    stest = "FC",
                    abs = TRUE,
                    kernel = TRUE,
                    n.perm = 1000,
                    correction =TRUE)
             standardGeneric("gseaK") )

#' @aliases gseaK-method
#' @rdname gseaK-methods
setMethod("gseaK", signature("data.frame"),
         function(expr, pheno, gSets, stest, abs, kernel, n.perm, correction){

  ####expression data
  expr <- data.matrix(expr)[,-1]
  gene_labels<-rownames(expr)
  sample_names<-colnames(expr)

  ###phenotype data
  class <- unlist(strsplit(pheno[[3]], " "))
  n_cl <-length(class)  #number of samples
  class<-as.factor(class)
  levels(class)=list('0'='c','1'='d')

  ####gene sets
  gSets<-data.frame(gSets[-1,])

  ## sample division into 2 groups
  n_gr<-as.numeric(table(class)) ### size of the groups
  in_gr1=which(class==0)
  in_gr2=which(class==1)

 ###### Rank metrics
 # if(stest=="FC")
    p_stat_t<-as.matrix(diffmean.stat(t(expr),class))

  rownames(p_stat_t)=rownames(expr)

  if (abs=="TRUE") {
    stat<-as.data.frame(abs(p_stat_t[,1]))
  } else {
    stat<-as.data.frame(p_stat_t[,1])}

  ###ranking of genes
  ord=stat[order(stat,decreasing = T),,drop=F]
  gene_name_ord<-vector(mode="character", length=0)
  gene_name_ord<-rownames(ord)


  n_perm<-1000 ###number of permutations
  n_gene_set<-ncol(gSets)
  result<-matrix(nrow = n_gene_set,ncol=8)
  fraction<-matrix(nrow=1,ncol=1)
  colnames(result)<-c('Gene Set','ES_obs','p-value','NES','FDR','','p_val_kerel','p-val_BH_ker')


  ###shuffle
  n_genes=nrow(expr)
  index<-seq(1,n_cl,1)

  index_ph<-permutation(index,n_perm)

t=1
  for (t in 1:n_gene_set){
    da=as.matrix(gSets[,t][!is.na(gSets[,t])])  ##Gene set
    name_gene<-as.vector(colnames(gSets[t]))  ##name of GS
    splitt<-strsplit(name_gene,'path.')
    name_gene_set<-splitt[[1]][2]
    colnames(da)<-name_gene_set


    N_H<-nrow(da)  #n in GS
    N<-length(gene_name_ord)  #N ordered genes

    #####Genes in GS
    l<-ginGS(gene_name_ord, as.character(da), ord[,1]);

    N_R<-sum(abs(as.numeric((l$stat))))
    P_miss<-(1/(N-N_H))

    ES_matrix<- ES(ord[,1], gene_name_ord, P_miss,N_R,l$poz)


    colnames(ES_matrix)=c("P_hit","P_miss","ES")


    rownames(ES_matrix)<-gene_name_ord

    sum_ES=max(abs(ES_matrix[,3]),na.rm = T)
    x_ES<-which(abs(ES_matrix[,3])==sum_ES)

    ES_obs<-ES_matrix[x_ES,3]

    ######...plots

    ES_p=matrix(nrow=n_perm)

    ###permutations

    for(i in 1:n_perm){
    expr2 <- expr[,index_ph[i,]]
    p_stat_t_p<-as.matrix(diffmean.stat(t(expr2),class))

    p_stat_t_p<-as.data.frame(abs(p_stat_t_p))


    ord_p=p_stat_t_p[order(p_stat_t_p,decreasing = T),,drop=F]   ###ranking of genes
    gene_name_ord_p<-as.matrix(rownames(ord_p))

    N<-dim(ord_p)[1]

    #####Genes in GS
    l_p<-ginGS(gene_name_ord_p, as.character(da), ord_p[,1]);

    N_R_p<-sum(abs(as.numeric((l_p$stat))))

    ES_matrix_p<-ES(ord_p[,1], rownames(ord_p), P_miss,N_R_p,l_p$poz)

    rownames(ES_matrix_p)<-gene_name_ord_p

   max_ES_p <- max(abs(ES_matrix_p[,3]),na.rm = T)
   x_ES_p<-which(abs(ES_matrix_p[,3])==max_ES_p)
   ES_p[i]<-ES_matrix_p[x_ES_p,3]  ## obserwowane ES z kazdej permtacji
    }


    ###plots


    ### p-value

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

    result[t,1]<-name_gene_set
    result[t,2]<-ES_obs
    result[t,3]<-p_val
#    result[t,4]<-NES ##NES observed
  ##  result[t,6]<-fraction
    result[t,7]<-p_val_kern
  }

  ##### Benjamini and Hochberg p-value correction


  ###poprawić!!!!
 # FDR <-function(x){
  #  for (p in 1:dim(x)[1]){
   #   if (x[p,3]>=0){
  #      NS<-length(which(x[,4]>=0))
  #      di<-length(which(x[-p,4]>=x[p,4]))
  #    }else{
  #      NS<-length(which(x[,4]<=0))
  #      di<-length(which(x[-p,4]<=x[p,4]))
  #    }
  #    x[p,5]=(as.numeric(x[p,6])*NS)/di
  #  }
  #  return(x)
#  }

 # FDR_results<-FDR(result)
#  FDR_results[,8]<-p.adjust(FDR_results[,7], method="BH")

 # res<-FDR_results[,-6]


  #### save results
 # write.csv(res,file = "wyniki_balanced_smyth.gmx",quote=F,row.names = F,sep="\t")
#  saveRDS(res,file='wynik_balanced_smyth.RData')

  new("matrix", result)

         })
