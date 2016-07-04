setwd("c:/docum180612/Francesc/fosses")
# csv file with haplotypes in rows and frequencies in columns. No row or column headers
freqalleles<-read.csv("xhaps.csv", header=FALSE)

# Set manually the number of loci and maximum number of different alleles/haplotypes
maxk<-4
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
sexcl<-0
nlr<-0
cexcl<-0
for (r in rr){
lr<-0
excl<-0
slr<-0
nexcl<-0

for (k in kk){

  vmare1<-rmultinom(1,1,freqalleles[,k])
  vmare2<-rmultinom(1,1,freqalleles[,k])
  vpare<-rmultinom(1,1,freqalleles[,k])
  
  for (kkk in aa){
    if (vmare1[kkk] == 1) {mare1<-kkk}
    if (vmare2[kkk]== 1){ mare2<-kkk}
    if (vpare[kkk]== 1){ pare<-kkk}
  
  }
  filla1<-pare
  dau <-runif(1)
  if (dau<0.5) {filla2<-mare1} else {filla2<-mare2}
  
  if (falsepaternity==1){
  filla1<-mare1
  filla2<-mare2
  }
  
  s<-0
  t<-0
  
  #dropout
  #if the paternal allele drops out, the whole locus is ignored in the calculations
  dau<-runif(1)
  if (dau>dp){
  s<-1
  excl<-0
  if (pare==filla1 & pare==filla2){t<-1/freqalleles[pare,k]} 
  if ((pare==filla1 & pare!=filla2) |(pare!=filla1 & pare==filla2)) {t<-1/(2*freqalleles[pare,k])}
  if (pare!=filla1 & pare!=filla2){
    excl<-1
  nexcl<-nexcl+1}
  
  if (s>0 & excl==0){slr<-slr+log10(t)}
  
}

if (s>0 & nexcl==0) {lr<-lr+slr}

slr<-0
}
sexcl<-sexcl+nexcl
if (s>0 & nexcl>0){cexcl<-cexcl+1}

if (nexcl==0){nlr<-nlr+1
vlr[nlr]<-lr}

}

if (nlr>0){
  vvlr<-matrix(0:0,nrow=nlr)
  vvlr<-vlr[1:nlr]}
if (nlr==0) {vvlr<-vlr}

mexcl<-sexcl/nrep

meanexcl<-sexcl/nrep
meanvlr<-mean(vvlr)
twopointfive <- quantile(vvlr, probs=0.025)
medianlr <-quantile(vvlr,probs=0.5)
ninesevenpointfive <- quantile (vvlr, probs=0.975)
print (paste0('chance of exclusion=',cexcl/nrep))
print (paste0('mean number of exclusionary loci=',mexcl))

#prints LR statistics for the subset of cases in which there was no exclusions

print (paste0('mean LR=',round(meanvlr,4)))
print (twopointfive)
print (medianlr)
print (ninesevenpointfive)





# Each single value of the log10(LR) is written out to a csv file. Edit the file name as necessary
write.csv(vlr,file="a_parefillaXfalsedp05.csv")