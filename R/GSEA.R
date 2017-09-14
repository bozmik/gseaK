install.packages('st')
library(st)

data("exprData")
data("phenoData")
data("geneSets")

expr <- exprData
pheno <- phenoData
gSets <- geneSets

gseaK <- function(expr, pheno, gSets, stest, abs, kernel, n.perm, correction) {

  ####expression data
  expr <- data.matrix(expr)[,-1]
  gene_labels<-rownames(expr)
  sample_names<-colnames(expr)

  ###cls file
  class <- unlist(strsplit(pheno[[3]], " "))
  n_cl <-length(class)  #number of samples
  class<-as.factor(class)
  levels(class)=list('0'='c','1'='d')

  ####gmx file - gene sets
  gSets<-data.frame(gSets[-1,])

  ## sample division into 2 groups
  n_gr<-as.numeric(table(class)) ### size of the groups
  in_gr1=which(class==0)
  in_gr2=which(class==1)

 ###### Rank metrics
  if(stest=="FC")
    p_stat_t<-as.matrix(diffmean.stat(t(expr),class))

  rownames(p_stat_t)=rownames(expr)


  if (abs=="TRUE") {
    stat<-as.data.frame(abs(p_stat_t[,1]))
  } else {
    stat<-as.data.frame(p_stat_t[,1])}

  ###ranking of genes
  ord=stat[order(stat,decreasing = T),,drop=F]
  gene_name_ord<-as.matrix(rownames(ord))


  n_perm<-1000 ###number of permutations
  n_gene_set<-ncol(gSets)
  result<-matrix(nrow = n_gene_set,ncol=8)
  fraction<-matrix(nrow=1,ncol=1)
  colnames(result)<-c('Gene Set','ES_obs','p-value','NES','FDR','','p_val_kerel','p-val_BH_ker')


  for (t in 1:n_gene_set){
    da=as.matrix(gSets[,t][!is.na(gSets[,t])])  ##Gene set
    name_gene<-as.vector(colnames(gSets[t]))  ##name of GS
    splitt<-strsplit(name_gene,'path.')
    name_gene_set<-splitt[[1]][2]
    colnames(da)<-name_gene_set


    N_H<-nrow(da)  #n in GS
    N<-length(gene_name_ord)  #N ordered genes

    #####Genes in GS

    names_GS<-matrix(nrow=N,ncol=1)
    matrix_stat<-matrix(nrow=N,ncol=1)
    set<-matrix(nrow=N,ncol=1)

    #####!!!!!!poprawić!!!
    for(k in 1:N){
      wh2<-which(as.vector(gene_name_ord[k]==da,'numeric')==1)
      if (length(wh2)!=0){
        matrix_stat[k,]<-ord[k,1]
        names_GS[k,1]<-gene_name_ord[k]
        set[k]=k
      }
    }

    in_GS<-cbind(names_GS,matrix_stat)
    in_GS<-na.omit(in_GS) #which genes are in gene set
    colnames(in_GS)<-c("gene name","statistic")
    set<-na.omit(set)

    N_R<-sum(abs(as.numeric((in_GS[,2]))))
    P_miss<-(1/(N-N_H))

    data_input<-list(ord=ord,gene_name_ord=gene_name_ord,P_miss=P_miss,N_R=N_R,N=N,da=da)

    ##fun_ES_compiled <- cmpfun(fun_ES)

    #ES_matrix<-fun_ES_compiled(data_input)

    ES_matrix<- ES(data_input)


    rownames(ES_matrix)<-gene_name_ord

    sum_ES=max(abs(ES_matrix[,5]),na.rm = T)
    x_ES<-which(abs(ES_matrix[,5])==sum_ES)

    ES_obs<-ES_matrix[x_ES,5]

    ######...plots

    n_ph=length(class)
    ES_p<-matrix(nrow=n_perm,ncol=1)
    NES_p<-matrix(nrow=n_perm,ncol=1)

    n_genes=nrow(A)
    index<-seq(1,n_ph,1)


    ###permutations

    ####poprawić!!!! w C++?

   ##cl <- makeCluster(detectCores(), type='PSOCK')
  ##  registerDoParallel(cl)

  ##ES_p<-foreach(i=1:n_perm,.combine=c) %dopar% perm_compiled(A,phenotyp,ph,P_miss,da=da,f_abs=f_abs)
  ##  ES_p<-unlist(ES_p)

  ##  stopCluster(cl)

    ###plots


    ### p-value

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


    NES_p<-fun(ES_p)
    NES_p<-cbind(ES_p,NES_p)
    colnames(NES_p)<-c('ES_p','NES_p')


    if(NES>=0){
      fr<-length(which(NES_p[,2]>=NES))  #fraction
      fr_poz<-length(which(NES_p[,2]>=0))
      fraction=fr/fr_poz
    }else{
      fr<-length(which(NES_p[,2]<=NES))
      fr_neg<-length(which(NES_p[,2]<0))
      fraction=fr/fr_neg
    }

    result[t,1]<-name_gene_set
    result[t,2]<-ES_obs
    result[t,3]<-p_val
    result[t,4]<-NES ##NES observed
    result[t,6]<-fraction
    result[t,7]<-p_val_kern

  }




  ##### Benjamini and Hochberg p-value correction


  ###poprawić!!!!
  FDR <-function(x){
    for (p in 1:dim(x)[1]){
      if (x[p,3]>=0){
        NS<-length(which(x[,4]>=0))
        di<-length(which(x[-p,4]>=x[p,4]))
      }else{
        NS<-length(which(x[,4]<=0))
        di<-length(which(x[-p,4]<=x[p,4]))
      }
      x[p,5]=(as.numeric(x[p,6])*NS)/di
    }
    return(x)
  }

  FDR_results<-FDR(result)
  FDR_results[,8]<-p.adjust(FDR_results[,7], method="BH")

  res<-FDR_results[,-6]


  #### save results
  write.csv(res,file = "wyniki_balanced_smyth.gmx",quote=F,row.names = F,sep="\t")
  saveRDS(res,file='wynik_balanced_smyth.RData')


}
