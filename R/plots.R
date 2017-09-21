# ---------------------------------------------------------------------------- #
#' plot Phit+Pmiss
#'
#' @param ES_matrix \code{matrix} that stores Phit, Pmiss, and ES distributions
#' @return Phit and Pmiss distributions.
#'
#' @export
######
plotDis<-function(ES_matrix){
plot(ES_matrix[,1],lwd=3,type = "l",col='blue',ylab =c('Cumulative distribution'),xlab = '')
lines(ES_matrix[,2],col='red',lwd=3)
legend('bottomright',c(expression('P'['miss']),expression('P'['hit'])),col=c('red','blue'),cex=1.2,pch=16,bty = "n")
}
# ---------------------------------------------------------------------------- #
#' plot ES
#'
#' @param ES_matrix \code{matrix} that stores ES distribution
#'
#' @return plot of ES distribution
#'
#' @export
######
plotES<-function(ES_matrix, ES_obs, x_ES, pos){
name_gene_set <-c('a1')
name_file<-paste(name_gene_set,c('.png'),sep="", collapse="")
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
