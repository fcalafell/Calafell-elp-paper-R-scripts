setwd("c:/docum180612/Francesc/fosses/def")
# csv file with alleles in rows and STRs in columns. No row or column headers
freqalleles<-read.csv("alleles.csv", header=FALSE)

# Set manually the number of loci and maximum number of different alleles
maxk<-121
maxa<-61

#Edit the number of iterations if needed
nrep=10000

#Edit the dropout fraction if needed 
dp<-0.5

#Set hetonly to 0 or 1 as needed
hetonly<-1

#Set falsepaternity to 0 or 1 as needed
falsepaternity<-1

kk<-c(1:maxk)
aa<-c(1:maxa)
rr<-c(1:nrep)
vlr<-matrix(0:0,nrow=nrep)

drawallele <- function (x, freqalleles){
  vecallele<-matrix(0:0,nrow=maxa)
  vecallele<-rmultinom(1,1,freqalleles[,x])
  
  for (kkk in aa){
    if (vecallele[kkk] == 1) {allele<-kkk}
     }
  return(allele)
}



for (r in rr){
lr<-1
for (k in kk){
  s<-999
  a1 <- drawallele(k,freqalleles)
  a2 <- drawallele(k,freqalleles)
  a3 <- drawallele(k,freqalleles)
  a4 <- drawallele(k,freqalleles)
  n2 <- drawallele(k,freqalleles)
  
 
 t1 <-a2
 t2 <- a4
 g1 <-a2
 g2 <- a4
 
  dau1 <- runif(1)
  if (dau1 < 0.5){t1<-a1} 
  
  dau2 <- runif(1)
  if (dau2 < 0.5){t2<-a3} 

  dau3 <- runif(1)
  if (dau3 < 0.5){g1<-a1} 
  
  dau4 <- runif(1)
  if (dau4 < 0.5){g2<-a3} 
  
 n1<- g2
  dau5 <- runif(1)
    if (dau5 < 0.5){n1<-g1}

if (falsepaternity==1){
t1<-a1
t2<-a2
n1<-a3
n2<-a4
}


 dau6<-runif(1)
 if (dau6<dp){t1<-t2}
 
  
 if (((t1 != n1) & (t1 != n2)) & ((t2 != n1) & (t2 != n2))){
   s<-0
 }

 
 if(((t1==t2)&(n1==n2))&(t1==n1)){
   s<-1/(2*freqalleles[t1,k])
 }

 if (((t1==t2) & (n1!=n2)) & ((t1==n1)|(t1==n2))) {
   s<-1/(4*freqalleles[t1,k])
 } 

 if (((n1==n2) & (t1!=t2)) & ((t1==n1)|(t2==n1))) {
   s<-1/(4*freqalleles[n1,k])
 } 

 if (((n1!=n2) & (t1!=t2)) & ((t1==n1)|(t2==n1))) {
   s<-1/(4*freqalleles[n1,k])
 }
 
 if (((t1!=t2) & (n1!=n2)) & ((t1==n1)|(t2==n1))){
   s<-1/(8*freqalleles[n1,k])
 }
 
 if (((t1!=t2) & (n1!=n2)) & ((t1==n2)|(t2==n2))){
   s<-1/(8*freqalleles[n2,k])
 }

 if (((t1!=t2) & (n1!=n2)) & (((t1==n1)&(t2==n2)) | ((t1==n2) & (t2==n1)))){
   s<-(1/(8*freqalleles[n1,k]))+(1/(8*freqalleles[n2,k]))
 }

if (hetonly=1){
 if (t1==t2){s<-0.5}
}

lr<-lr*(0.5+s)

}
vlr[r]<-lr
}


meanvlr<-mean(vlr)
twopointfive <- quantile(vlr, probs=0.025)
medianlr <-quantile(vlr,probs=0.5)
ninesevenpointfive <- quantile (vlr, probs=0.975)
# The mean, median, and 95% CI of LR are printed out

print (meanvlr)
print (twopointfive)
print (medianlr)
print (ninesevenpointfive)

# Each single value of the LR is written out to a csv file. Edit the file name as necessary
write.csv(vlr,file="avinetfalsedp05.csv")