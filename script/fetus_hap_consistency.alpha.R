args=commandArgs(T)
sList=args[1]
workdir=args[2]
tag="alpha"

# HBA
gene.chr="chr16"
gene.start=c(222846)
gene.end=c(227520)
gene.name="HBA1-HBA2"

extL=2e5
xl=c(0,4e5)

col.unmatch="grey"
col.matchF="lightseagreen"
col.matchM="lightcoral"

#pdf(paste0(workdir,"/","FigureS2_",tag,".pdf"),width=16,height=25)
tiff(paste0(workdir,"/","FigureS2_",tag,".600dpi.tiff"),width=16,height=25,units='in',res=600,compression='jpeg',family='Arial')
par(mar=c(5,2,0,0),mgp=c(3,.5,0),mfrow=c(1,2))
read.table(sList,stringsAsFactors = F)[,1]->flist
nS=length(flist)

#mat
plot(0,t='n',
     xlab="Position (chromosome 16)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=xl,ylim=c(0.5,nS+1.7),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(0,xl[2],by=1e5)
#axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)
xlab=round(seq(0,3.5e5,by=5e4)/1e6,digits=2)
xlab=sprintf("%.2f", xlab)
axis(1,at=seq(0,3.5e5,by=5e4),labels=xlab,lwd=0.75,cex.axis=1.3)
for(i in 1:nS){
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.mout.plot"),stringsAsFactors = T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchM
    t[t$pos<=gene.end[1]+extL&t$pos>=0,]->t
    t$start=0
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i-0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nS-i-0.15,nS-i-0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    points(t$pos,rep(nS-i+0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nS-i+0.15,nS-i+0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    text(3e4,nS-i,labels =paste0(flist[i]),pos=2,cex=1.5)
}
rect(gene.start,-0.5,gene.end,nS-0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nS-0.15,labels = gene.name,cex=1.5)
text(1.9e5,nS+0.5,labels ="Maternal haplotype",cex=2,font=2)

#pat
plot(0,t='n',
     xlab="Position (chromosome 16)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=xl,ylim=c(0.5,nS+1.7),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(0,xl[2],by=1e5)
#axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)
xlab=round(seq(0,3.5e5,by=5e4)/1e6,digits=2)
xlab=sprintf("%.2f", xlab)
axis(1,at=seq(0,3.5e5,by=5e4),labels=xlab,lwd=0.75,cex.axis=1.3)
for(i in 1:nS){
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.fout.plot"),stringsAsFactors = T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchF
    t[t$pos<=gene.end[1]+extL&t$pos>=0,]->t
    t$start=0
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i-0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nS-i-0.15,nS-i-0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    points(t$pos,rep(nS-i+0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nS-i+0.15,nS-i+0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    text(3e4,nS-i,labels =paste0(flist[i]),pos=2,cex=1.5)
}
rect(gene.start,-0.5,gene.end,nS-0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nS-0.15,labels = gene.name,cex=1.5)
text(1.9e5,nS+0.5,labels ="Paternal haplotype",cex=2,font=2)

legend(1.8e5,nS+2.8,legend=c("Paternal matched SNPs","Maternal matched SNPs","Parental mismatched SNPs"),pch=20,cex=1.5,col=c(col.matchF,col.matchM,col.unmatch))
dev.off()
