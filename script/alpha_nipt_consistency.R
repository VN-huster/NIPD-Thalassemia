args=commandArgs(T)
tiff(paste0("FigureS3_alpha.600dpi",".tiff"),width=20,height=35,units='in',res=600,compression='jpeg',family='Arial')
gene=c(222846,227520)
gene.name="HBA1-HBA2"
extL=2e5
par(family="Times")
par(mar=c(8,6,5,1),mgp=c(3,.5,0),mfrow=c(1,2))
xl=c(0,3.5e5)
y=30

read.table("alpha_relation.list.sort",stringsAsFactors = F)->flist
plot(0,t='n',
     xlim=xl,ylim=c(10,(nrow(flist)+1)*y),
     axes=F,
     ylab="",
     xlab="Position (chromosome 16)/Mbp",
     cex.lab=1.5,
     font=2,
     yaxs='i')
title_size=1.4
for(i in nrow(flist):1){
    rect(xl[1],-10+i*y,xl[2],10+i*y,col='gray93',density=-1,border=F)
    lb=seq(-10,10,by=10)
    xlab=round(seq(0,3.5e5,by=5e4)/1e6,digits=2)
    xlab=sprintf("%.2f", xlab)
    xlab <- paste0(xlab,"")
    ylab=seq(-10,10,by=10)
    axis(2,at=ylab+i*y,labels=c(-10,0,10),xaxp=c(-9+i*y,10+i*y,18),las=1,cex.axis=1.2)

    polygon(x=c(xl[1],xl[2]),y=c(0+i*y,0+i*y),density=-1,border='grey',lwd=2) 
    points(x=c(1,1)*xl[1],y=c(5,-5)+i*y,cex=c(4,4)*1,col=c('black','black'),pch=c(21,21),bg=c('white','white'))
    text(x=xl[1],y=c(5,-5)+i*y,labels=c("P","N"),col=c('black','black'),cex=1.5,font=2)
    if (file.exists(paste0(flist[i,1],"_nipt.OR.mout"))){
	read.table(paste0(flist[i,1],"_nipt.OR.mout"),header=T,stringsAsFactors=F)->a
	a->tmp2
	if (!is.null(tmp2$OddRatio_M)){
	    points(tmp2$pos,tmp2$OddRatio_M+i*y,xlim=xl,t='l',pch=20,col="lightcoral",cex=1,lwd=3,lty=1)
	    points(tmp2$pos,tmp2$OddRatio_M+i*y,xlim=xl,t='p',pch=20,col="lightcoral",cex=1,lwd=1,lty=1)
	    mtext(paste0("m",flist[i,1]),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y,)
	}else{
	    if (tmp2$P0>0.5){
		points(tmp2$pos,tmp2$P0/tmp2$P1+i*y,xlim=xl,t='p',pch=20,col="lightcoral",cex=1,lwd=3,lty=1)
	    }else{
		points(tmp2$pos,-tmp2$P0/tmp2$P1+i*y,xlim=xl,t='p',pch=20,col="lightcoral",cex=1,lwd=3,lty=1)
	    }
	    mtext(paste0("m",flist[i,1],"(#)"),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y,)
	}
    }else{
	mtext(paste0("m",flist[i,1],"(#)"),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y,)
    }
}
axis(1,at=seq(0,3.5e5,by=5e4),labels=xlab,lwd=0.75,cex.axis=title_size)
mtext("Maternal inheritance",side=3,col=1,cex=2,line=0,font=2)
mtext(expression(paste('Odds ratio in plasma / ',log[e],' (P/N)')),side=2, line=3.8,col=1,outer=F,cex=1.5,at=18*y,)
rect(gene[1],0,gene[2],nrow(flist)*y+10,col= NA, border='dark violet',lty=3,lwd=2)
text(mean(gene[1],gene[2]),nrow(flist)*y+20,labels = gene.name,cex=title_size,font=2)

plot(0,t='n',
     xlim=xl,ylim=c(10,(nrow(flist)+1)*y),
     axes=F,
     ylab="",
     xlab="Position (chromosome 16)/Mbp",
     cex.lab=1.5,
     font=2,
     yaxs='i')
title_size=1.4
for(i in nrow(flist):1){
    rect(xl[1],-10+i*y,xl[2],10+i*y,col='gray93',density=-1,border=F)
    lb=seq(-10,10,by=10)
    xlab=round(seq(0,3.5e5,by=5e4)/1e6,2)
    xlab=sprintf("%.2f", xlab)
    xlab <- paste0(xlab,"")
    ylab=seq(-10,10,by=10)
    axis(2,at=ylab+i*y,labels=c(-10,0,10),xaxp=c(-9+i*y,10+i*y,18),las=1,cex.axis=1.2)
    polygon(x=c(xl[1],xl[2]),y=c(0+i*y,0+i*y),density=-1,border='grey',lwd=2)
    points(x=c(1,1)*xl[1],y=c(5,-5)+i*y,cex=c(4,4)*1,col=c('black','black'),pch=c(21,21),bg=c('white','white'))
    text(x=xl[1],y=c(5,-5)+i*y,labels=c("P","N"),col=c('black','black'),cex=1.5,font=2)
    if (file.exists(paste0(flist[i,1],"_nipt.OR.fout"))){
	read.table(paste0(flist[i,1],"_nipt.OR.fout"),header=T,stringsAsFactors=F)->a
	a->tmp2
	if (!is.null(tmp2$OddRatio_F)){
	    if (nrow(tmp2)==1 && tmp2$OddRatio_F==0){
		if (tmp2$P0/tmp2$P1<10){
		    if (tmp2$P0>0.5){
    			points(tmp2$pos,tmp2$P0/tmp2$P1+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
		    }else{
			points(tmp2$pos,-tmp2$P0/tmp2$P1+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
		    }
		}else{
		    if (tmp2$P0>0.5){
    			points(tmp2$pos,10+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
		    }else{
			points(tmp2$pos,-10+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
		    }
		}
		mtext(paste0("p",flist[i,1],"(#)"),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
	    }else{
    		points(tmp2$pos,tmp2$OddRatio_F+i*y,xlim=xl,t='l',pch=20,col="lightseagreen",cex=1,lwd=3,lty=1)
		points(tmp2$pos,tmp2$OddRatio_F+i*y,xlim=xl,t='p',pch=20,col="lightseagreen",cex=1,lwd=1,lty=1)
		mtext(paste0("p",flist[i,1]),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
	    }
	}
    }else{
	mtext(paste0("p",flist[i,1],"(#)"),side=2, line=2.3,col=1,outer=F,cex=1.5,at=i*y)
    }
}
axis(1,at=seq(0,3.5e5,by=5e4),labels=xlab,lwd=0.75,cex.axis=title_size)
mtext("Paternal inheritance",side=3,col=1,cex=2,line=0,font=2)
mtext(expression(paste('Odds ratio in plasma / ',log[e],' (P/N)')),side=2, line=3.8,col=1,outer=F,cex=1.5,at=18*y,)
rect(gene[1],0,gene[2],nrow(flist)*y+10,col= NA, border='dark violet',lty=3,lwd=2)
text(mean(gene[1],gene[2]),nrow(flist)*y+20,labels = gene.name,cex=title_size,font=2)

dev.off()




