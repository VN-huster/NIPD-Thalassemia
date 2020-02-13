args=commandArgs(T)
sList=args[1]
workdir=args[2]
tag="beta"

# HBA
gene.chr="chr11"
gene.start=c(5246696)
gene.end=c(5248301)
gene.name="HBB"

extL=5e5
xl=c(gene.start-extL,5.5e6)

col.unmatch="grey"
col.matchF="lightseagreen"
col.matchM="lightcoral"

pdf(paste0(workdir,"/","FigureS2_",tag,".pdf"),width=16,height=25)
par(mar=c(5,2,0,0),mgp=c(3,.5,0),mfrow=c(1,2))
read.table(sList,stringsAsFactors = F)[,1]->flist
nS=length(flist)

#mat
plot(0,t='n',
     xlab="Position (chromosome 11)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=c(4.7e6,5.5e6),ylim=c(0.5,nS+1.7),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(4.7e6,5.5e6,by=1e5)
axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)
for(i in 1:nS){
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.mout.plot"),stringsAsFactors = T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchM
    t[t$pos<=gene.end[1]+extL&t$pos>=gene.start[1]-extL,]->t
    t$start=gene.start[1]-extL
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i-0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+2e4,xl[2]-1e5),y=c(nS-i-0.15,nS-i-0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    points(t$pos,rep(nS-i+0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+2e4,xl[2]-1e5),y=c(nS-i+0.15,nS-i+0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    text(xl[1],nS-i,labels =paste0(flist[i]),pos=2,cex=1.5)
}
rect(gene.start,-0.5,gene.end,nS-0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nS-0.15,labels = gene.name,cex=1.5)
text(5.1e6,nS+0.5,labels ="Maternal haplotype",cex=2,font=2)

#pat
plot(0,t='n',
     xlab="Position (chromosome 11)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=c(4.7e6,5.5e6),ylim=c(0.5,nS+1.7),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(4.7e6,5.5e6,by=1e5)
axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)

for(i in 1:nS){
    read.table(paste0(workdir,"/",flist[i],"_verify.",tag,".OUT.stat.fout.plot"),stringsAsFactors = T,header=T)->t
    t$match=col.unmatch
    nrow(t)
    t[t$parent==t$fetus,]$match=col.matchF
    t[t$pos<=gene.end[1]+extL&t$pos>=gene.start[1]-extL,]->t
    t$start=gene.start[1]-extL
    t$end=gene.end[1]+extL
    points(t$pos,rep(nS-i-0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+2e4,xl[2]-1e5),y=c(nS-i-0.15,nS-i-0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    points(t$pos,rep(nS-i+0.15,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
    polygon(x=c(xl[1]+2e4,xl[2]-1e5),y=c(nS-i+0.15,nS-i+0.15),col='grey65',density=-1,border='grey65',lwd=.5)
    text(xl[1],nS-i,labels =paste0(flist[i]),pos=2,cex=1.5)
}
rect(gene.start,-0.5,gene.end,nS-0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nS-0.15,labels = gene.name,cex=1.5)
text(5.1e6,nS+0.5,labels ="Paternal haplotype",cex=2,font=2)

legend(xl[2]-4e5,nS+2.5,legend=c("Paternal matched SNPs","Maternal matched SNPs","Parental mismatched SNPs"),pch=20,cex=1.5,col=c(col.matchF,col.matchM,col.unmatch))
dev.off()
