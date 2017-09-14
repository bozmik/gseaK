
fun <- function(x){
  NES_k<-matrix(nrow=length(x),ncol=1)
  for (l in 1:length(x)){
    if (x[l]>=0){
      m<-mean(x[x>=0])
      NES_k[l]=x[l]/m}
    else{
      m<-mean(abs(x[x<0]))
      NES_k[l]=x[l]/m
    }
  }
  return(NES_k)
}

#fun_compiled <- cmpfun(fun)

fun <- function(x){
  NES_k<-matrix(nrow=length(x),ncol=1)
  for (l in 1:length(x)){
    if (x[l]>=0){
      m<-mean(x[x>=0])
      NES_k[l]=x[l]/m}
    else{
      m<-mean(abs(x[x<0]))
      NES_k[l]=x[l]/m
    }
  }
  return(NES_k)
}

#fun_compiled <- cmpfun(fun)
