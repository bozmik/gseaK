#### in C++!!!!
ES <- function(x){
  P_hit_vec<-matrix(nrow=x$N)
  P_miss_vec<-matrix(nrow=x$N)
  ES<-matrix(nrow=x$N)
  is_in_GS<-matrix(nrow=x$N)
  ES_matrix<-cbind(x$ord[,1],is_in_GS,P_hit_vec,P_miss_vec,ES)
  colnames(ES_matrix)=c("t-statistic,","presence in GS","P_hit","P_miss","ES")

  for (k in 1:N){
    wh<-x$gene_name_ord[k]==x$da
    wh2<-which(wh=="TRUE")
    if (length(wh2)!=0){
      if (k==1){
        ES_matrix[k,2]=1
        sum_t_stat=abs(ES_matrix[k,1])
        P_hit<-sum_t_stat/x$N_R
        ES_matrix[k,3]=P_hit
        ES_matrix[k,4]=0
      }else{
        ES_matrix[k,2]=1
        sum_t_stat=abs(ES_matrix[k,1])
        P_hit<-sum_t_stat/x$N_R
        ES_matrix[k,3]=ES_matrix[k-1,3]+P_hit
        ES_matrix[k,4]=ES_matrix[k-1,4]
      }
    }else if(length(wh2)==0){
      if (k==1) {
        ES_matrix[k,2]=0
        ES_matrix[k,3]=0
        ES_matrix[k,4]=x$P_miss
      }else{
        ES_matrix[k,2]=0
        ES_matrix[k,3]=ES_matrix[k-1,3]
        ES_matrix[k,4]=ES_matrix[k-1,4]+x$P_miss
      }
    }
    ES_matrix[k,5]=ES_matrix[k,3]-ES_matrix[k,4]
  }
  return(ES_matrix)
}





