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
#' @param plot1 plot1
#' @param plot2 plot2
#' @param plot3 plo3
#' @return returns a \code{matrix} object
#' @useDynLib gseaK
#' @importFrom Rcpp sourceCpp
#' @docType methods
#' @examples
#'data("exprData")
#'data("phenoData")
#'data("geneSets")
#'
#'gseaK(expr = exprData, pheno = phenoData, gSets = geneSets[,1:10])
#'gseaK(expr = exprData, pheno = phenoData, gSets = geneSets[,1:10], stest = "shrinkage.t")
#'
#' @rdname gseaK-methods
#' @export
setGeneric("gseaK",
           function(expr, pheno, gSets,
                    stest = "FC",
                    abs = TRUE,
                    kernel = TRUE,
                    n.perm = 1000,
                    correction =TRUE,
                    plot1 = TRUE,
                    plot2 = TRUE,
                    plot3 = TRUE)
             standardGeneric("gseaK") )

#' @aliases gseaK-method
#' @rdname gseaK-methods
setMethod("gseaK", signature("data.frame"),
         function(expr, pheno, gSets, stest, abs, kernel, n.perm, correction, plot1, plot2, plot3){

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
  gSets[1,]<-colnames(gSets)

  ## sample division into 2 groups
  n_gr<-as.numeric(table(class)) ### size of the groups
  in_gr1=which(class==0)
  in_gr2=which(class==1)

 ###### Rank metrics
  if(stest == "FC"){
    p_stat_t <- as.matrix(diffmean.stat(t(expr),class))
  } else if(stest == "shrinkage.t"){
   p_stat_t <- as.matrix(shrinkt.stat(t(expr),class,var.equal=F,paired=FALSE, verbose=T))
  } else if(stest == "2s.Bayesian"){
    p_stat_t <- as.matrix(RankingFoxDimmic(expr,class,type = "unpaired")@statistic)
  } else if(stest == "Efron"){
    p_stat_t <- as.matrix(RankingEbam(expr,class,type = "unpaired")@statistic)
  } else if(stest == "SAM"){
    p_stat_t <- as.matrix(RankingSam(expr,class,type = "unpaired")@statistic)
  } else if(stest == "penalized.t"){
    p_stat_t <- as.matrix(RankingSoftthresholdT(expr,class, type = c("unpaired"))@statistic)
  } else if(stest == "bayesian.t"){
     p_stat_t <- as.matrix(RankingBaldiLong(expr,class, type = c("unpaired"))@statistic)
  } else if(stest == "moderatet.t"){
    p_stat_t <- as.matrix(RankingLimma(expr,class,type = c("unpaired"))@statistic)
  }

  ##
  ##else if(stest == "moderated.wt"){
  ##  p_stat_t <- as.matrix(mwt(expr,class,log.it = FALSE)$MWT)
  ##}

  rownames(p_stat_t)=rownames(expr)

  if (abs==TRUE) {
    stat<-as.data.frame(abs(p_stat_t[,1]))
  } else {
    stat<-as.data.frame(p_stat_t[,1])}

  ###ranking of genes
  ord=stat[order(stat,decreasing = T),,drop=F]
  gene_name_ord<-vector(mode="character", length=0)
  gene_name_ord<-rownames(ord)


  n_gene_set<-ncol(gSets)
  result<-matrix(nrow = 7,ncol=n_gene_set)
  fraction<-matrix(nrow=1,ncol=1)


  ###shuffle
  n_genes=nrow(expr)
  index<-seq(1,n_cl,1)

  index_ph<-permutation(index,n.perm)

  result <- sapply(gSets, GSet, gene_name_ord, ord[,1], n.perm, expr, class, index_ph, stest, abs, plot1, plot2, plot3, simplify=TRUE)
  rownames(result)<-c('ES_obs','p-value','NES','FDR','','p_val_kerel','p-val_BH_ker')

  ##### Benjamini and Hochberg p-value correction

  ##result[4,] <- FDR(result[5,], result[3,])

  result[7,]<-p.adjust(result[6,], method="BH")

 # res<-FDR_results[,-6]


  #### save results
 # write.csv(res,file = "wyniki_balanced_smyth.gmx",quote=F,row.names = F,sep="\t")
#  saveRDS(res,file='wynik_balanced_smyth.RData')

  new("matrix", result)

         })
