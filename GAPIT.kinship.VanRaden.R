`GAPIT.kinship.VanRaden` <-
function(snps,hasInbred=TRUE) {
# Object: To calculate the kinship matrix using the method of VanRaden (2009, J. Dairy Sci. 91:4414???C4423)
# Authors: Zhwiu Zhang
# Last update: october 25, 2014 
############################################################################################## 

print("Calculating kinship with VanRaden method...")

#snps=hm$GD 
nSNP=ncol(snps)
nInd=nrow(snps)
n=nInd 
snpMean= apply(snps,2,mean)   #get mean for each snp
print("substracting mean...")
snps=t(snps)-snpMean    #operation on matrix and vector goes in direction of column
print("Getting X'X...")
#K=tcrossprod((snps), (snps))
K=crossprod((snps), (snps)) #Thanks to Peng Zheng, Meng Huang and Jiafa Chen for finding the problem
print("Adjusting...")
#Extract diagonals
i =1:n
j=(i-1)*n
index=i+j
d=K[index]
DL=min(d)
DU=max(d)
floor=min(K)

K=(K-floor)/(DL-floor)
MD=(DU-floor)/(DL-floor)

#Handler of diagonals over 2
#print("MD")
#print(MD)
#print(K[1:5,1:5])

if(is.na(K[1,1])) stop ("GAPIT says: Missing data is not allowed for numerical genotype data")
if(MD>2)K[index]=K[index]/(MD-1)+1

#Handler of inbred
if(MD<2 & hasInbred) K=2*K/((DU-floor)/(DL-floor))

print("Calculating kinship with VanRaden method: done")

return(K)
}

