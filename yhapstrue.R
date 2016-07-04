setwd("c:/docum180612/Francesc/fosses/def")
# csv input file with haplotypes as rows, last column is haplotype absolute frequency

inputtable<-read.csv("yhaps.csv", header=FALSE)
# Set manually the number of loci, number of haplotypes, and number of invididuals
nloci<-14
nhap<-867
nind<-1003
niter<-10000

#Edit the dropout fraction if needed
dp<-0.5

nn<-c(1:nhap)
ii<-c(1:niter)
jj<-c(1:nloci)

lr<-matrix(0:0,nrow=niter)
haphap<-matrix(0:0, nrow=nloci)
freqhap<-matrix(0:0,nrow=nhap)
freqhap<-inputtable[,nloci+1]
freqhaprel<-freqhap/nind
slr<-0
for (i in ii){
  dau<-rmultinom(1,1,freqhaprel)
  for (n in nn){
    if (dau[n]==1){hap<-n}
  }
  nmiss<-0
  for (j in jj){
    haphap[j]=inputtable[hap,j]
    dau2<-runif(1)
    if (dau2<dp){
      haphap[j]<- -99
  nmiss<-nmiss+1}
  }
  freq<-0
  
  if (nmiss==0) {lr[i]<- (-1)*log10(freqhaprel[hap])} else {
  for (n in nn){
    matx<-0
    for (j in jj){
      if (haphap[j]==inputtable[n,j]){matx<-matx+1}
      if (haphap[j]==-99){matx<-matx+1}
    }
    if (matx==nloci){freq<-freq+freqhaprel[n]}
  }
  
  lr[i]<- (-1)*log10(freq)}
  
}
mlr<-mean(lr)
# Mean log10(LR) is printed
print (mlr)

# log10(LR) values are written to a file. Change name as necessary.
write.csv(lr,file="ytruedp05.csv")