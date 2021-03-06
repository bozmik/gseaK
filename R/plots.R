# ---------------------------------------------------------------------------- #
#' plot Phit+Pmiss
#'
#' @param ES_matrix \code{matrix} that stores Phit, Pmiss, and ES distributions
#' @param name_gene_set GS name
#' @return Phit and Pmiss distributions.
#'
#' @export
######
plotDis<-function(ES_matrix, name_gene_set){
name_file<-paste(name_gene_set,c('plot1'),c('.png'),sep="_", collapse = "")
png(filename = name_file)

plot(ES_matrix[,1],lwd=3,type = "l",col='blue',ylab =c('Cumulative distribution'),xlab = '', main=name_gene_set)
lines(ES_matrix[,2],col='red',lwd=3)
legend('bottomright',c(expression('P'['miss']),expression('P'['hit'])),col=c('red','blue'),cex=1.2,pch=16,bty = "n")
dev.off()
}

# ---------------------------------------------------------------------------- #
#' plot ES
#'
#' @param ES_matrix \code{matrix} that stores ES distribution
#' @param ES_obs ES value
#' @param x_ES position of ES value
#' @param pos positions of genes from GS
#' @param name_gene_set GS name
#' @return plot of ES distribution
#'
#' @export
######
plotES<-function(ES_matrix, ES_obs, x_ES, pos, name_gene_set){
name_file<-paste(name_gene_set,c('plot2'),c('.png'),sep="_", collapse = "")
png(filename = name_file)
image_gene <- par(mfrow=c(2, 1))
N <- dim(ES_matrix)[1]
par(mfrow=c(2,1))
par(mar=c(0, 4.5, 1, 0.5))  #bottom, left, top, and right.
plot(1:N,as.vector(ES_matrix[,3]),axes=F,pch=20,ylab='Enrichment Score',xlab=NA,ylim =c(min(ES_matrix[,3])-0.1,1),type="l",col='blue',lwd = 3,main=name_gene_set)
lines(c(1:N),c(rep(0,N)),lwd=1,lty=1,col="black")
lines(c(x_ES,x_ES),c(0,ES_obs),col='red',lty=10,lwd=2)
text(x_ES,ES_obs+0.1, round(ES_obs, 2),cex=1)
Axis(side=2,las = 1, tck = -0.02)
par(mar=c(5.5,4.5,2, 0.5))
plot(pos,rep(0.1,length(pos)),axes=F,type='h',col='blue',ylim=c(0,0.1),xlim=c(1,N),ylab =NA,xlab=NA)
lines(rep(0,N),lwd=1.5,lty=1,col="black")
Axis(side=1,las = 1, tck = -0.1)
mtext(side = 1, "Gene List Rank", line = 2.5)
#savePlot(filename=name_file,type="png")
dev.off()
}

# ---------------------------------------------------------------------------- #
#' plot of ES_null
#'
#' @param ES_p \code{vector} that stores ES_null values
#' @param ES_obs ES value
#' @param dens Kernel density estimation
#' @param name_gene_set GS name
#' @return plot of ES_null
#'
#' @export

plotPerm<-function(ES_p, ES_obs, dens, name_gene_set){
if (min(ES_p)>ES_obs){
  x_lim1=ES_obs-0.1
  x_lim2=max(ES_p)+0.1
} else if(max(ES_p)<ES_obs) {
  x_lim2=ES_obs+0.1
  x_lim1=min(ES_p)-0.1
}else {
  x_lim1=min(ES_p)-0.1
  x_lim2=max(ES_p)+0.1}

name_file<-paste(name_gene_set,c('plot3'),c('.png'),sep="_", collapse = "")
png(filename = name_file)
r<-hist(ES_p,freq=F,breaks=100,main=name_gene_set,xlim=c(x_lim1,x_lim2),xlab='ES')
rug(ES_p)
#lines(ES_p)
lines(dens,col='red',lwd = 2)
lines(c(ES_obs, ES_obs), c(0, max(dens$y)-0.2), col = "steelblue", lwd = 3, lty = 22)
text(ES_obs,(max(dens$y)-0.1), round(ES_obs, 4),cex=1.5)
legend("topright",c(c("kernel density estymation"),c("observed ES")),pch=16,col=c("red","steelblue"),cex=1,bty = "n")
dev.off()
}
