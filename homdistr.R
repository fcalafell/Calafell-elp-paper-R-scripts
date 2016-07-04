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

kk<-c(1:maxk)
aa<-c(1:maxa)
rr<-c(1:nrep)
hom<-matrix(0:0,nrow=nrep)

for (r in rr){
shom<-0
for (k in kk){

  val1<-rmultinom(1,1,freqalleles[,k])
  val2<-rmultinom(1,1,freqalleles[,k])
  
  for (kkk in aa){
    if (val1[kkk] == 1) {al1<-kkk}
    if (val2[kkk]== 1){ al2<-kkk}
    }
  dau<-runif(1)
  if (dau<dp){al1<-al2}
    if (al1==al2){shom<-shom+1}
  
}
hom[r]<-shom
}

meanhom<-mean(hom)
twopointfive <- quantile(hom, probs=0.025)
medianhom <-quantile(hom,probs=0.5)
ninesevenpointfive <- quantile (hom, probs=0.975)

# prints the mean, median, and 95% CI of expected number of homozygotes

print (meanhom)
print (twopointfive)
print (medianhom)
print (ninesevenpointfive)

#writes the whole distribution of the number of expected homozygotes. Please rename file as needed
write.csv(hom,file="homdp05.csv")