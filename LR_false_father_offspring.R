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

shet<-0
sexcl<-0
cexcl<-0
nlr<-0
sumlr<-0
for (r in rr){
lr<-1
excl<-0
nexcl<-0
nhet<-0
for (k in kk){

  va1<-rmultinom(1,1,freqalleles[,k])
  va2<-rmultinom(1,1,freqalleles[,k])
  vb1<-rmultinom(1,1,freqalleles[,k])
  vb1<-rmultinom(1,1,freqalleles[,k])
  
  for (kkk in aa){
    if (va1[kkk] == 1) {a1<-kkk}
    if (va2[kkk]== 1){ a2<-kkk}
    if (vb1[kkk] == 1){ b1<-kkk}
    if (vb1[kkk] == 1){ b2<-kkk}
      }
 
  dau <- runif(1)
  if (dau < dp) {a2 <- a1}
 
if (hetonly==0 | (hetonly== 1 & a1!=a2)){
  nhet<-nhet+1
  if ((a1!=b1 & a2!=b2) | (a2!=b1 & a1!=b2)){
    excl<-1
    nexcl<-nexcl+1
  }
  else {
    if (excl==0){
    if (a1==b1){
      ac<-a1
      afnc<-b2
    apnc<-a2}
    
    if (a1==b2){
      ac<-a1
      afnc<-b1
      apnc<-a2}
    
    if (a2==b1){
      ac<-a2
      afnc<-b2
      apnc<-a1}
    
    if (a2==b2){
      ac<-a2
      afnc<-b1
      apnc<-a1}
    



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
}


shet<-shet+nhet
sexcl<-sexcl+nexcl
if (excl==0) {nlr<-nlr+1
              vlr[nlr]<-lr
              sumlr<-sumlr+lr} else {cexcl<-cexcl+1}
}



if (nlr>0){
  vvlr<-matrix(0:0,nrow=nlr)
  vvlr<-vlr[1:nlr]}
if (nlr==0) {vvlr<-vlr}

meanhet<-shet/nrep
meanexcl<-sexcl/nrep
meanvlr<-mean(vvlr)
twopointfive <- quantile(vvlr, probs=0.025)
medianlr <-quantile(vvlr,probs=0.5)
ninesevenpointfive <- quantile (vvlr, probs=0.975)
print (paste0('chance of exclusion=',cexcl/nrep))

#prints LR statistics for the subset of cases in which there was no exclusions
print (paste0('mean LR=',round(sumlr/nlr,4)))
print (twopointfive)
print (medianlr)
print (ninesevenpointfive)
write.csv(vvlr,file="test.csv")