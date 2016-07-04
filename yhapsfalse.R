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
haphap1<-matrix(0:0, nrow=nloci)
haphap2<-matrix(0:0, nrow=nloci)
hapcons<-matrix(0:0, nrow=nloci)
freqhap<-matrix(0:0,nrow=nhap)
freqhap<-inputtable[,nloci+1]
freqhaprel<-freqhap/nind

cexcl<-0
clr<-0
for (i in ii){
  dau<-rmultinom(1,1,freqhaprel)
  for (n in nn){
    if (dau[n]==1){hap1<-n}
  }
  dau<-rmultinom(1,1,freqhaprel)
  for (n in nn){
    if (dau[n]==1){hap2<-n}
  }
  
  nmiss1<-0
  for (j in jj){
    haphap1[j]=inputtable[hap1,j]
    dau2<-runif(1)
    if (dau2<dp){
      haphap1[j]<- -99
      nmiss1<-nmiss1+1}
  }
  nmiss2<-0
  for (j in jj){
    haphap2[j]=inputtable[hap2,j]
    dau2<-runif(1)
    if (dau2<dp){
      haphap2[j]<- -99
      nmiss2<-nmiss2+1}
  }

  if (nmiss1+nmiss2==0 & hap1==hap2){
    clr<-clr+1
    lr[clr]<-log10(freqhap[hap1])
  }
  
  if (nmiss1+nmiss2==0 & hap1!=hap2){
   cexcl<-cexcl+1
  }
  
    ndif<-0
    for (j in jj){
      if ((haphap1[j]!=haphap2[j]) & (haphap1[j]!=-99) & (haphap2[j]!=-99)) {ndif=ndif+1}
    }
    if (ndif>0) {cexcl<-cexcl+1} else {
      for (j in jj){
        if ((haphap1[j]==-99) & (haphap2[j]==-99)) {hapcons[j]<- -99}
        if ((haphap1[j]==-99) & (haphap2[j]!=-99)) {hapcons[j]<- haphap2[j]}
        if ((haphap1[j]!=-99) & (haphap2[j]==-99)) {hapcons[j]<- haphap1[j]}
        if ((haphap1[j]!=-99) & (haphap2[j]!=-99)) {hapcons[j]<- haphap1[j]}
      }
      freq<-0
      for (n in nn){
        matx<-0
        for (j in jj){
          if (hapcons[j]==inputtable[n,j]){matx<-matx+1}
          if (hapcons[j]==-99){matx<-matx+1}
        }
        if (matx==nloci){freq<-freq+freqhaprel[n]}
    }
    if (freq==0){freq<-1/(nind+1)}
    clr<-clr+1
    lr[clr]<- (-1)*log10(freq)
  }
}


if (clr>0){
  vvlr<-matrix(0:0,nrow=clr)
  vvlr<-lr[1:clr]}
if (clr==0) {vvlr<-lr}


meanvlr<-mean(vvlr)
twopointfive <- quantile(vvlr, probs=0.025)
medianlr <-quantile(vvlr,probs=0.5)
ninesevenpointfive <- quantile (vvlr, probs=0.975)
print (paste0('chance of exclusion=',cexcl/niter))

#prints log10(LR) statistics for the subset of cases in which there was no exclusions
print (paste0('mean LR=',round(meanvlr,4)))
print (twopointfive)
print (medianlr)
print (ninesevenpointfive)

# log10(LR) values are written to a file. Change name as necessary.
write.csv(vvlr,file="test.csv")