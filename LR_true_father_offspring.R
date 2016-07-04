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

kk<-c(1:maxk)
aa<-c(1:maxa)
rr<-c(1:nrep)
vlr<-matrix(0:0,nrow=nrep)
dp<-0.5
for (r in rr){
lr<-1
for (k in kk){

  vac<-rmultinom(1,1,freqalleles[,k])
  vapnc=rmultinom(1,1,freqalleles[,k])
  vafnc=rmultinom(1,1,freqalleles[,k])
  
  for (kkk in aa){
    if (vac[kkk] == 1) {ac<-kkk}
    if (vapnc[kkk]== 1){ apnc<-kkk}
    if (vafnc[kkk] == 1){ afnc<-kkk}
  }
 
  dau <- runif(1)
  if (dau < dp) apnc <- ac
if (hetonly=0){
 if (ac!=apnc) {

  if ((ac==apnc) & (ac==afnc)) {
    t<-1
    s<-0}
  
  if ((ac==apnc) & (ac!=afnc)) {
    t<-2
    s<-0}
  
  if ((ac!=apnc) & (ac==afnc)) {
    t<-2
    s<-0}
  
  if ((ac!=apnc) & (ac!=afnc)  & (apnc!=afnc)) {
   
      t<-4
      s<-0
    }
    
  if ((ac!=apnc) & (ac!=afnc)  & (apnc==afnc)){
    t<-4
    s<-1/(4*freqalleles[apnc,k])
  }
  lr<-lr*((1/(t*freqalleles[ac,k]))+s)  
  }
}
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
write.csv(vlr,file="test.csv")