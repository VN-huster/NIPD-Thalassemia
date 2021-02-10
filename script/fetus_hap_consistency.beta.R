args=commandArgs(T)
sList=args[1]
workdir=args[2]
tag="beta"

gene.chr="chr11"
gene.start=c(5246696)
gene.end=c(5248301)
gene.name="HBB"

extL=5e5
xl=c(4.7e6,5.4e6)

col.unmatch="grey"
col.matchF="lightseagreen"
col.matchM="lightcoral"

tiff(paste0(workdir,"/","FigureS2_",tag,".600dpi.tiff"),width=13,height=18,units='in',res=600,compression='jpeg',family='Arial')
par(las=1,mar=c(5,4,0,0))
read.table(sList,stringsAsFactors = F)[,1]->flist
nS=length(flist)

plot(0,t='n',
     xlab="Position (chromosome 11)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=c(4.7e6,5.5e6),ylim=c(0.5,nS+1.7),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(4.7e6,5.5e6,by=1e5)
xlab=round(seq(xl[1],xl[2],by=(xl[2]-xl[1])/7)/1e6,2)
xlab=sprintf("%.2f", xlab)
axis(1,at=seq(xl[1],xl[2],by=(xl[2]-xl[1])/7),labels=xlab,lwd=0.75,cex.axis=1.3)

for(i in 1:nS){
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.fout.plot"),stringsAsFactors = T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchF
    t[t$pos<=gene.end[1]+extL&t$pos>=gene.start[1]-extL,]->t
    t$start=gene.start[1]-extL
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i-0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+4e4,xl[2]),y=c(nS-i-0.15,nS-i-0.15),border='grey65',density=-1,lwd=1)
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.mout.plot"),stringsAsFactors =T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchM
    t[t$pos<=gene.end[1]+extL&t$pos>=gene.start[1]-extL,]->t
    t$start=gene.start[1]-extL
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i+0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+4e4,xl[2]),y=c(nS-i+0.15,nS-i+0.15),border='grey65',density=-1,lwd=1)
    text(xl[1]+3e4,nS-i,labels =paste0(flist[i]),pos=2,cex=1.5)
}
rect(gene.start,-0.5,gene.end,nS-0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nS-0.15,labels = gene.name,cex=1.5)

legend(xl[2]-2.5e5,nS+2.5,legend=c("Maternal matched SNPs","Paternal matched SNPs","Parental mismatched SNPs"),pch=20,cex=1.5,col=c(col.matchM,col.matchF,col.unmatch))
dev.off()
