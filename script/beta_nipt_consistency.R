args=commandArgs(T)
tiff(paste0("FigureS3_beta.600dpi",".tiff"),width=20,height=30,units='in',res=600,compression='jpeg',family='Arial')
gene=c(5246696,5248301)
gene.name="HBB"
extL=2.5e5
par(family="Times")
par(mar=c(8,6,5,1),mgp=c(3,.5,0),mfrow=c(1,2))
xl=c(4.7e6,5.4e6)
y=30

read.table("beta_relation.list.sort",stringsAsFactors = F)->flist
plot(0,t='n',
     xlim=xl,ylim=c(10,(nrow(flist)+1)*y),
     axes=F,
     ylab="",
     xlab="Position (chromosome 11)/Mbp",
     cex.lab=1.5,
     font=2,
     yaxs='i')
title_size=1.4
for(i in nrow(flist):1){
    read.table(paste0(flist[i,1],"_nipt.OR.mout"),header=T,stringsAsFactors=F)->a
    a->tmp2
    rect(xl[1],-10+i*y,xl[2],10+i*y,col='gray93',density=-1,border=F)
    points(tmp2$pos,tmp2$OddRatio_M+i*y,xlim=xl,t='l',pch=20,col="lightcoral",cex=1,lwd=3,lty=1);
    points(tmp2$pos,tmp2$OddRatio_M+i*y,xlim=xl,t='p',pch=20,col="lightcoral",cex=1,lwd=1,lty=1);
    lb=seq(-10,10,by=10)
    xlab=round(seq(xl[1],xl[2],by=(xl[2]-xl[1])/7)/1e6,2)
    xlab=sprintf("%.2f", xlab)
    xlab <- paste0(xlab,"")
    ylab=seq(-10,10,by=10)
    axis(2,at=ylab+i*y,labels=c(-10,0,10),xaxp=c(-9+i*y,10+i*y,18),las=1,cex.axis=1.2)
    mtext(paste0("m",flist[i,1]),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
    polygon(x=c(xl[1],xl[2]),y=c(0+i*y,0+i*y),density=-1,border='grey',lwd=2)
    points(x=c(1,1)*xl[1],y=c(5,-5)+i*y,cex=c(4,4)*1,col=c('black','black'),pch=c(21,21),bg=c('white','white'));
    text(x=xl[1],y=c(5,-5)+i*y,labels=c("P","N"),col=c('black','black'),cex=1.5,font=2)
}
axis(1,at=seq(xl[1],xl[2],by=(xl[2]-xl[1])/7),labels=xlab,lwd=0.75,cex.axis=title_size)
mtext("Maternal inheritance",side=3,col=1,cex=2,line=0,font=2)
mtext(expression(paste('Odds ratio in plasma / ',log[e],' (P/N)')),side=2, line=3.8,col=1,outer=F,cex=1.5,at=13*y,)
rect(gene[1],0,gene[2],nrow(flist)*y+10,col= NA, border='dark violet',lty=3,lwd=2)
text(mean(gene[1],gene[2]),nrow(flist)*y+20,labels = gene.name,cex=title_size,font=2)

plot(0,t='n',
     xlim=xl,ylim=c(10,(nrow(flist)+1)*y),
     axes=F,
     ylab="",
     xlab="Position (chromosome 11)/Mbp",
     cex.lab=1.5,
     yaxs='i')
title_size=1.4
for(i in nrow(flist):1){
    read.table(paste0(flist[i,1],"_nipt.OR.fout"),header=T,stringsAsFactors=F)->a
    a->tmp2
    rect(xl[1],-10+i*y,xl[2],10+i*y,col='gray93',density=-1,border=F)
    points(tmp2$pos,tmp2$OddRatio_F+i*y,xlim=xl,t='l',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
    points(tmp2$pos,tmp2$OddRatio_F+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=1,lty=1)
    lb=seq(-10,10,by=10)
    xlab=round(seq(xl[1],xl[2],by=(xl[2]-xl[1])/7)/1e6,2)
    xlab=sprintf("%.2f", xlab)
    xlab <- paste0(xlab,"")
    ylab=seq(-10,10,by=10)
    axis(2,at=ylab+i*y,labels=c(-10,0,10),xaxp=c(-9+i*y,10+i*y,18),las=1,cex.axis=1.2)
    if (flist[i,1]=="F23"){
	mtext(paste0("p",flist[i,1],"(#)"),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
    }else{
	mtext(paste0("p",flist[i,1]),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
    }
    polygon(x=c(xl[1],xl[2]),y=c(0+i*y,0+i*y),density=-1,border='grey',lwd=2)
    points(x=c(1,1)*xl[1],y=c(5,-5)+i*y,cex=c(4,4)*1,col=c('black','black'),pch=c(21,21),bg=c('white','white'));
    text(x=xl[1],y=c(5,-5)+i*y,labels=c("P","N"),col=c('black','black'),cex=1.5,font=2)
}
axis(1,at=seq(xl[1],xl[2],by=(xl[2]-xl[1])/7),labels=xlab,lwd=0.75,cex.axis=title_size)
mtext("Paternal inheritance",side=3,col=1,cex=2,line=0,font=2)
mtext(expression(paste('Odds ratio in plasma / ',log[e],' (P/N)')),side=2, line=3.8,col=1,outer=F,cex=1.5,at=13*y,)
rect(gene[1],0,gene[2],nrow(flist)*y+10,col= NA, border='dark violet',lty=3,lwd=2)
text(mean(gene[1],gene[2]),nrow(flist)*y+20,labels = gene.name,cex=title_size,font=2)

dev.off()
