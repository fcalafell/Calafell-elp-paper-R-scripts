setwd("c:/docum180612/Francesc/fosses")
# csv file with haplotypes in rows and frequencies in columns. No row or column headers
freqalleles<-read.csv("xhaps.csv", header=FALSE)

# Set manually the number of loci and maximum number of different loci
maxk<-4
# Set manually the number of loci and maximum number of different alleles/haplotypes per locus
maxa<-106

#Edit the number of iterations if needed
nrep=10000

#Edit the dropout fraction if needed 
dp<-0.5

#Set falsepaternity to 0 or 1 as needed
falsepaternity<-1


kk<-c(1:maxk)
aa<-c(1:maxa)
rr<-c(1:nrep)
vlr<-matrix(0:0,nrow=nrep)
for (r in rr){
lr<-0
for (k in kk){

  vmare1<-rmultinom(1,1,freqalleles[,k])
  vavi<-rmultinom(1,1,freqalleles[,k])
  
  
  for (kkk in aa){
    if (vmare1[kkk] == 1) {mare1<-kkk}
    if (vavi[kkk]== 1){ avi<-kkk}
    
  
  }
  mare2<-avi
  dau <-runif(1)
  if (dau<0.5) {net<-mare1} else {net<-mare2}
  
  if (falsepaternity==1){net<-mare1}
  
  s<-0
  t<-0
  #dropout
  #if the paternal allele drops out, the whole locus is ignored in the calculations
  dau<-runif(1)
  if (dau>dp){
  s<-0.5  
  if (net==avi){t<-1/(2*freqalleles[net,k])} else{t<-0}
  
  
  
  if (s>0){lr<-lr+log10(s+t)}  
}
}



vlr[r]<-lr
}
meanvlr<-mean(vlr)
twopointfive <- quantile(vlr, probs=0.025)
medianlr <-quantile(vlr,probs=0.5)
ninesevenpointfive <- quantile (vlr, probs=0.975)
# The mean, median, and 95% CI of the log10(LR) are printed out
print (meanvlr)
print (medianlr)
print (twopointfive)
print (ninesevenpointfive)


# Each single value of the log10(LR) is written out to a csv file. Edit the file name as necessary
write.csv(vlr,file="O_netmaternfalsedp05.csv")