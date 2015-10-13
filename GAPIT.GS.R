`GAPIT.GS` <-
function(KW,KO,KWO,GAU,UW){
#Object: to derive BLUP for the individuals without phenotype
#Output: BLUP
#Authors: Zhiwu Zhang 
# Last update: April 17, 2011 
##############################################################################################
#print(length(UW))
UO=try(t(KWO)%*%solve(KW)%*%UW)
if(inherits(UO, "try-error")) UO=t(KWO)%*%ginv(KW)%*%UW

n=ncol(UW) #get number of columns, add additional for individual name

#Assign BLUP of group to its individuals
BLUP=data.frame(as.matrix(GAU[,1:4]))
#BLUP.W=BLUP[which(GAU[,3]==1),]
BLUP.W=BLUP[which(GAU[,3]<2),]
order.W=order(as.numeric(as.matrix(BLUP.W[,4])))
ID.W=as.numeric(as.matrix(BLUP.W[order.W,4]))
n.W=max(ID.W)
DS.W=diag(n.W)[ID.W,]
ind.W=DS.W%*%UW
all.W=cbind(BLUP.W[order.W,],ind.W)
all=all.W

BLUP.O=BLUP[which(GAU[,3]==2),]
if(nrow(BLUP.O)>0){
order.O=order(as.numeric(as.matrix(BLUP.O[,4])))
ID.O=as.numeric(as.matrix(BLUP.O[order.O,4]))
n.O=max(ID.O)
DS.O=diag(n.O)[ID.O,]
ind.O=DS.O%*%UO
all.O=cbind(BLUP.O[order.O,],ind.O)
all=rbind(all.W,all.O)
}

colnames(all)=c("Taxa", "Group", "RefInf","ID","BLUP","PEV")

#print("GAPIT.GS accomplished successfully!")
return(list(BLUP=all))
}#The function GAPIT.GS ends here

