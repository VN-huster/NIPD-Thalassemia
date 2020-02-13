args=commandArgs(T)
prefix=args[1]


initHMM<-function (States, Symbols, startProbs = NULL, transProbs = NULL,
    emissionProbs = NULL)
{
    nStates = length(States)
    nSymbols = length(Symbols)
    S = rep(1/nStates, nStates)
    T = 0.5 * diag(nStates) + array(0.5/(nStates), c(nStates,
        nStates))
    E = array(1/(nSymbols), c(nStates, nSymbols))
    names(S) = States
    dimnames(E) = list(states = States, symbols = Symbols)
    if (!is.null(startProbs)) {
        S[] = startProbs[]
    }
    if (!is.null(transProbs)) {
        if(length(dim(transProbs))==2){
		T = array(0,c(dim(transProbs),1))
		T[,,1] = transProbs
	        dimnames(T) = list(from = States, to = States,NULL)
        }else if(length(dim(transProbs))==3){
                T = transProbs
	        dimnames(T) = list(from = States, to = States,NULL)
        }
    }
    if (!is.null(emissionProbs)) {
        E[, ] = emissionProbs[, ]
    }
    return(list(States = States, Symbols = Symbols, startProbs = S,
        transProbs = T, emissionProbs = E))
}

viterbi<-function (hmm, observation)
{
    hmm$transProbs[is.na(hmm$transProbs)] = 0
    hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    nObservations = length(observation)
    nStates = dim(hmm$transProbs)[1]
    if(length(dim(hmm$transProbs))==3 &&dim(hmm$transProbs)[3]!=nObservations-1){
        cat("Warning:The 3rd dim of transProbs(",dim(hmm$transProbs)[3],")
                is not as long as observations-1(",nObservations-1,")\n")
    }

    tritransProbs = array(rep(c(hmm$transProbs),round(1+nStates**2*nObservations/length(hmm$transProbs))),
                                c(nStates,nStates,nObservations-1))
    dimnames(tritransProbs) = dimnames(hmm$transProbs)
    nStates = length(hmm$States)
    v = array(NA, c(nStates, nObservations))
    dimnames(v) = list(states = hmm$States, index = 1:nObservations)
    for (state in hmm$States) {
        v[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,
            observation[1]])
    }
    if(nObservations>=2){
    for (k in 2:nObservations) {
        for (state in hmm$States) {
            maxi = NULL
            for (previousState in hmm$States) {
                temp = v[previousState, k - 1] + log(tritransProbs[previousState,
                  state,k-1])
                maxi = max(maxi, temp)
            }
            v[state, k] = log(hmm$emissionProbs[state, observation[k]]) +
                maxi
        }
    }}
    optimal_path_prob.temp<- max(v[, nObservations])
    print (paste0("optimal_path_prob:",optimal_path_prob.temp))
    viterbiPath = rep(NA, nObservations)
    for (state in hmm$States) {
        if (max(v[, nObservations]) == v[state, nObservations]) {
            viterbiPath[nObservations] = state
            break
        }
    }
    if(nObservations>=2){
    for (k in (nObservations - 1):1) {
        for (state in hmm$States) {
            if (max(v[, k] + log(tritransProbs[, viterbiPath[k +
                1],k])) == v[state, k] + log(tritransProbs[state,
                viterbiPath[k + 1],k])) {
                viterbiPath[k] = state
                break
            }
        }
    }}
    return(viterbiPath)
}

# fetal fraction
read.table(paste0(prefix,".fra"),header=T,stringsAsFactors=F)->fra
as.numeric(fra$fra)->fra$fra
#fra[fra$fra>0,]->fra
boxplot.stats(fra$fra)->fra.box
fra[!fra$fra%in%fra.box$out,]->fra.fix
mean(fra.fix$fra)->fra.mean
sd(fra.fix$fra)->fra.sd
print(paste0("fetal fraction ",fra.mean,"(sd:",fra.sd,")"))
write.table(fra.mean,paste0(prefix,".concentration.txt"),sep="\t",quote=F,row.names=F,col.names=F)

t=fra$fra
t=t[t<0.5&t>0.01]
boxplot.stats(t)->t.box
t.fra=t.box$stats[3]
print(paste0("fetal fraction ",t.fra))
#fra.mean=t.fra



read.table(paste0(prefix,".OUT"),header=T)->OUT
OUT[OUT$chr!=24,]->OUT
OUT[OUT$F0==OUT$F1,]->md
OUT[OUT$M0==OUT$M1&OUT$F0!=OUT$F1,]->fd
fd[fd$chr!=23,]->fd
if(nrow(fd)==0){
	print ('No qualified snps')
	q()
}
#md[t.fr>t.fr.box$stats[2]&t.fr<t.fr.box$stats[4],]->md.fix
#md->md.fix

#md.fix[t.fr>t.fr.box$stats[2]&t.fr<t.fr.box$stats[4],]->md.fix



fd->fd.fix
uplimit=0.95
lowlimit=0.05

########################
#HMM
#######################
#genetic_map.txt="/vol6/home/bgi_wangyy/project/pl/genetic_map/genetic_map_GRCh37_all.txt";
genetic_map.txt="/THL4/home/bgi_guofengyu/work/haplotyping/config/genetic_map_GRCh37_HBB_HBA_1M.txt";
genetic_map <- read.table(genetic_map.txt,head=T)
getTransProbs  <- function(map,snpset,chrname){
	chrmap0<- map[map[,1]==chrname,c(2,4)]
	chrmap1<- cbind(snpset[snpset[,1]==chrname,2],-1)
	if(dim(chrmap0)[1]==0){
		if(chrname==23){chrname="X"}
		chrname <- paste("chr",chrname,sep="")
		chrmap0<- map[map[,1]==chrname,c(2,4)]
	}
	chrmap0<-as.matrix(chrmap0)
	if(length(chrmap0)<=2 || length(chrmap1)<=2){
		cat("chr ",chrname,length(chrmap0)/2,length(chrmap1)/2,"\n")
		NULL
	}else{
	chrmap0<- chrmap0[order(chrmap0[,1]),]
	chrmap1<- chrmap1[order(chrmap1[,1]),]
	n0     <- nrow(chrmap0)
	n1     <- nrow(chrmap1)
	i0     <- 1
	i1     <- 1
	cat("chr ",chrname,n0,n1,"\n")
	while(i1<=n1 && chrmap1[i1,1]<=chrmap0[1,1]){
		chrmap1[i1,2] <- 0
		i1 <- i1+1
	}
	while(i1<=n1){
		while(i0<=n0&&chrmap0[i0,1]<chrmap1[i1,1]){
			i0 <- i0+1
		}
		if(i0<=n0){
			chrmap1[i1,2] <- (chrmap0[i0,2]-chrmap0[i0-1,2])*(chrmap1[i1,1]-chrmap0[i0-1,1])/(chrmap0[i0,1]-chrmap0[i0-1,1])+chrmap0[i0-1,2]
			i1 <- i1+1
		}else{
			while(i1<=n1){
				chrmap1[i1,2] <- chrmap0[n0,2]
				i1 <- i1+1
			}
		}
	}
	TransProbs <- array(0,c(2,2,n1-1))
	for(i in 1:(n1-1)){
		p_ii1 <- (chrmap1[i+1,2]-chrmap1[i,2])/100
		if(p_ii1>0.5){p_ii1<-0.5}
		TransProbs[,,i] <- matrix(c(1-p_ii1, p_ii1, p_ii1, 1-p_ii1),2)
	}
	TransProbs
	}
}

#count=0
#sum=0
#for(i in 1:nrow(fd)){
#    t=fd[i,c("M_ref","M_alt")]
#    count=count+min(t)
#    sum=sum+sum(t)
#}
#Perror=count/sum
Perror=0.02

for(i in 1:nrow(fd.fix)){
    Mp_dep=fd.fix$Mp_ref[i]+fd.fix$Mp_alt[i]
    if(fd.fix$M0[i]==1){
        dbinom(fd.fix$Mp_ref[i],Mp_dep,Perror)->fd.fix$Hom[i]
        dbinom(fd.fix$Mp_ref[i],Mp_dep,0.5*fra.mean)->fd.fix$Het[i]
	if(fd.fix$F0[i]==1){ # Hom
	    fd.fix$Hom[i]/(fd.fix$Hom[i]+fd.fix$Het[i])->fd.fix$P0[i]
	}else{ # Het
	    fd.fix$Het[i]/(fd.fix$Hom[i]+fd.fix$Het[i])->fd.fix$P0[i]
	}
    }else if(fd.fix$M0[i]==0){
        dbinom(fd.fix$Mp_alt[i],Mp_dep,Perror)->fd.fix$Hom[i]
        dbinom(fd.fix$Mp_alt[i],Mp_dep,0.5*fra.mean)->fd.fix$Het[i]
	if(fd.fix$F0[i]==0){ # Hom
	    fd.fix$Hom[i]/(fd.fix$Hom[i]+fd.fix$Het[i])->fd.fix$P0[i]
	}else{ # Het
	    fd.fix$Het[i]/(fd.fix$Hom[i]+fd.fix$Het[i])->fd.fix$P0[i]
	}
    }

    fd.fix$P0[i]=min(fd.fix$P0[i],uplimit)
    fd.fix$P0[i]=max(fd.fix$P0[i],lowlimit)
    1-fd.fix$P0[i]->fd.fix$P1[i]
}

States=c("1","2")
startProbs=c(.5,.5)
#transProbs=matrix(c(.99,.01,.01,.99),2)

#paste0("chr",c(1:22,"X","Y"))->chr
c(1:22)->chr
upper=10
lower=-10
fd.fix$OddRatio_F<-'0'
for(i in 1:length(chr)){
    if(nrow(fd.fix[fd.fix$chr==chr[i],])>1){
        fd.fix[fd.fix$chr==chr[i],]->t
        1:nrow(t)->Symbols
        t(t[,c("P0","P1")])->emissionProbs
        getTransProbs(genetic_map,t,chr[i])->transProbs
        initHMM(States,Symbols,startProbs,transProbs,emissionProbs)->hmm
        1:nrow(t)->Observation
        viterbi(hmm,Observation)->t$fetal

		log(t$P0[1]/t$P1[1],10)->t$or[1]
		j=1
		transProbs[as.numeric(t$fetal[j]),1,j]->reP0
		transProbs[as.numeric(t$fetal[j]),2,j]->reP1
		log((t$P0[j]*reP0)/(t$P1[j]*reP1),10)->t$or[j]
		print (t$or[1])
		for(j in 1:(nrow(t)-1)){
			transProbs[as.numeric(t$fetal[j]),1,j]->reP0
			transProbs[as.numeric(t$fetal[j]),2,j]->reP1
			log((t$P0[j+1]*reP0)/(t$P1[j+1]*reP1),10)->t$or[j+1]
			t$or[j+1]=min(upper,max(lower,t$or[j+1]))
		}

        t$fetal->fd.fix[fd.fix$chr==chr[i],"fetal_F"]
		t$or->fd.fix[fd.fix$chr==chr[i],"OddRatio_F"]
    }
}
#("0")->fd.fix$fetal_M



write.table(fd.fix,paste0(prefix,".OR.fout"),sep="\t",quote=F,row.names=F)
