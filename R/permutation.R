perm<-function(A,phenotyp,ph,P_miss,da,f_abs){
  n_genes=nrow(A)

  n_gr<-as.numeric(table(class))


  n_ph=sum(n_gr)

  index<-seq(1,n_ph,1)
  index_ph<-sample(index,size = n_cl, replace = F)  ## random indexes for classes

  subset<- class[index_ph]

  #####statistical tests
  #####for #Welch two sample t-test,#Signal to Noise Ratio and #Wilcoxon rank sum test
  #in_gr1_p<-index_ph[1:n_gr[1]]
  #in_gr2_p<-index_ph[(n_gr[1]+1):(n_gr[1]+n_gr[2])]
  #p_stat_t_p=as.matrix(apply(A,1, function(x) test_compiled(x,index1=in_gr1_p,index2=in_gr2_p)))

  ####### for the rest of the methods
  A2=A[,index_ph]

  ######################## for Shrinkage t-statistic and FC
  #library(st)

  #### FC
  #p_stat_t_p<-as.matrix(diffmean.stat(t(A2),class))

  ####Shrinkage t-statistic
  # p_stat_t_p<-as.matrix(shrinkt.stat(t(A2),class,var.equal=F,paired=FALSE, verbose=T))

  #####################################
  #library(GeneSelector)

  ##### two-sample Bayesian t test
  #p_stat_t_p<-as.matrix(RankingFoxDimmic(A2,class,type = "unpaired")@statistic)

  ######Efron
  #p_stat_t_p<-as.matrix(RankingEbam(A2,class,type = "unpaired")@statistic)

  ####SAM
  #p_stat_t_p<-as.matrix(RankingSam(A2,class,type = "unpaired")@statistic)

  #####Penalized t-statistic
  #p_stat_t_p<-as.matrix(RankingSoftthresholdT(A2,class, type = c("unpaired"))@statistic)

  #####Bayesian t-statistic
  #p_stat_t_p<-as.matrix(RankingBaldiLong(A2,class, type = c("unpaired"))@statistic)

  #####moderatet t-test
  #p_stat_t_p<-as.matrix(RankingLimma(A2,class,type = c("unpaired"))@statistic)

  #####moderated Welch Test
  library(mwt)
  p_stat_t_p<-as.matrix(mwt(A2,class,log.it = FALSE)$MWT)

  rownames(p_stat_t_p)=row.names(A)


  if (f_abs==1) {
    p_stat_t_p<-as.data.frame(abs(p_stat_t_p))
  } else {
    p_stat_t_p<-as.data.frame(p_stat_t_p)}

  ord_p=p_stat_t_p[order(p_stat_t_p,decreasing = T),,drop=F]   ###ranking of genes
  gene_name_ord_p<-as.matrix(rownames(ord_p))

  N<-dim(ord_p)[1]

  names_GS_p<-matrix(nrow=N,ncol=1)
  matrix_stat_p<-matrix(nrow=N,ncol=1)
  set_p<-matrix(nrow=N,ncol=1)

  for(k in 1:N){
    wh2<-which(as.vector(gene_name_ord[k]==da,'numeric')==1)
    if (length(wh2)!=0){
      matrix_stat_p[k,]<-ord_p[k,1]
      names_GS_p[k,1]<-gene_name_ord_p[k]
      set_p[k]=k
    }
  }

  in_GS_p<-cbind(names_GS_p,matrix_stat_p)
  in_GS_p<-na.omit(in_GS_p) #which genes are in gene set
  colnames(in_GS_p)<-c("gene name","test")
  set_p<-na.omit(set_p)

  N_R_p<-sum(abs(as.numeric((in_GS_p[,2]))))


  P_hit_vec<-matrix(nrow=N)
  P_miss_vec<-matrix(nrow=N)
  ES<-matrix(nrow=N)
  is_in_GS<-matrix(nrow=N)
  ES_matrix_p<-cbind(ord_p,is_in_GS,P_hit_vec,P_miss_vec,ES)

  colnames(ES_matrix_p)=c("statistic,","presence in GS","P_hit","P_miss","ES")

  data_input_p<-list(ord=ord_p,gene_name_ord=rownames(ord_p),P_miss=P_miss,N_R=N_R_p,N=N,da=da)

  ES_matrix_p<-fun_ES_compiled(data_input_p)
  rownames(ES_matrix_p)<-gene_name_ord_p

  max_ES_p=max(abs(ES_matrix_p[,5]),na.rm = T)
  x_ES_p<-which(abs(ES_matrix_p[,5])==max_ES_p)
  ES_p<-ES_matrix_p[x_ES_p,5]  ## obserwowane ES z kazdej permtacji

  return(ES_p)
}



#perm_compiled<-cmpfun(perm)
