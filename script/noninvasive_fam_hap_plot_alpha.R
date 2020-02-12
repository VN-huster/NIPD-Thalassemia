args<-commandArgs(T)
read.table(args[1],header=T,stringsAsFactors=F)->a
mathet<-a[a$M0!=a$M1,]
pathet<-a[a$P0!=a$P1,]
meq<-0
mne<-0
peq<-0
pne<-0
peq<-0
pne<-0
for (i in 1:nrow(mathet)){
    if(mathet$M0[i]==mathet$Mathap1[i]){
	mathet$tag[i] <- 0.5 
	meq<-meq+1
    }else{
	mathet$tag[i] <- -0.5
	mne<-mne+1
    }
}
for (i in 1:nrow(pathet)){
    if(pathet$P0[i]==pathet$Pathap1[i]){
	pathet$tag[i] <- 0.5 
	peq<-peq+1
    }else{
	pathet$tag[i] <- -0.5
	pne<-pne+1
    }
}
pdf(paste0(args[2],"_famcomparison.pdf"),width=10,height=6)
par(mar=c(3,3,3,3), mfrow=c(2,1))
HBB_region<-c(222846, 227520)
xl<-c(0, HBB_region[2]+1e6)
yl<-c(-0.6,0.6)
plot(mathet$pos,mathet$tag, xlim=xl, ylim=yl, t='b', lwd=0.5, pch=20, col=2, xlab="Position (chr11)", ylab="",axes=F, main="Mat haplotype comparison")
axl<-seq(0,1.3e6,by=1e5)
xlab=round(axl/1e6,2)
xlab <- paste0(xlab,"Mb")
axis(1,at=axl,labels=xlab)
text(x=xl[2]-1e4,y=c(0.4,-0.4),labels=c(paste(meq,"SNPs"), paste(mne,"SNPs")),col='blue')
abline(v=HBB_region,col='light green')
#abline(v=c(215400,223300,219817),col=1)
abline(v=215400,col=1)
text(x=mean(HBB_region),y=0.6,labels="HBA2-HBA1")
abline(h=0,col=5)
plot(pathet$pos,pathet$tag, xlim=xl, ylim=yl, t='b', lwd=0.5, pch=20, col=2, xlab="Position (chr11)", ylab="",axes=F, main="Pat haplotype comparison")
axl<-seq(0,1.3e6,by=1e5)
xlab=round(axl/1e6,2)
xlab <- paste0(xlab,"Mb")
axis(1,at=axl,labels=xlab)
text(x=xl[2]-1e4,y=c(0.4,-0.4),labels=c(paste(peq,"SNPs"), paste(pne,"SNPs")),col='blue')
abline(v=HBB_region,col='light green')
#abline(v=c(215400,223300,219817),col=1)
abline(v=215400,col=1)
text(x=mean(HBB_region),y=0.6,labels="HBA2-HBA1")
abline(h=0,col=5)
