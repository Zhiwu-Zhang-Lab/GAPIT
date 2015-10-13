`GAPIT.Phenotype.Simulation` <-
function(GD,GM=NULL,h2=.75,NQTN=10,QTNDist="normal",effectunit=1,category=0){
#Object: To simulate phenotype from genotye
#Input: GD - n by m +1 dataframe or n by m big.matrix
#intput: h2 - heritability
#intput: NQTN - number of QTNs
#intput: QTNDist - Distribution of QTN, options are  "geometry", "normal"
#intput: effectunit - effect of fitst QTN, the nect effect is its squre
#intput: theSeed - seed for randomization
#Output: Y,U,E,QTN.Position, and effect
#Straitegy: NA
#Authors: Qishan Wang and Zhiwu Zhang
#Start  date: April 4, 2013
#Last update: April 4, 2013

#Set orientation
#Strategy: the number of rows in GD and GM are the same if GD has SNP as row

nm=ncol(GD)-1   #Initial by assume GD has snp in col
if(!is.null(GM)) nm=nrow(GM)
ngd1=nrow(GD)
ngd2=ncol(GD)
ngd1=abs(ngd1-nm)
ngd2=abs(ngd2-nm)
orientation="row"
ns=ncol(GD)
if(min(ngd1,ngd2)>0){
  orientation="col"
  ns=nrow(GD)
}



n= ns   #number of samples
m=nm  #number of markers

#Set QTN effects
if (QTNDist=="normal"){ addeffect<-rnorm(NQTN,0,1)
	}else
	{addeffect=effectunit^(1:NQTN)}


#Simulating Genetic effect
#r=sample(2:m,NQTN,replace=F)
QTN.position=sample(1:m,NQTN,replace=F)
if(orientation=="col") SNPQ=as.matrix(GD[,(QTN.position+1)])
if(orientation=="row") SNPQ=t(as.matrix(GD[QTN.position,]))

#Replace non-variant QTNs  (does not work yet)
#inComplete=TRUE
#while(inComplete){
#  inComplete=FALSE
#  myVar=apply(SNPQ,2,var)
#  index=which(myVar==0)
#  nInVar=length(index)
#  if(nInVar>0){
#    inComplete=TRUE
#    New.position=sample(1:m,nInVar,replace=F)
#    if(orientation=="col") SNPQ[,index]=as.matrix(GD[,(New.position+1)])
#    if(orientation=="row") SNPQ[,index]=t(as.matrix(GD[New.position,]))
#  }
#}#end of while


effect=SNPQ%*%addeffect
effectvar=var(effect)
residualvar=(effectvar-h2*effectvar)/h2

#Variance explained by each SNP
effectInd=SNPQ%*%diag(addeffect)
varInd=apply(effectInd,2,var)
effectSeq=order(varInd,decreasing = TRUE)

#Simulating Residual and phenotype
residual=rnorm(n,0,sqrt(residualvar))
if(orientation=="col") myY=cbind(as.data.frame(GD[,1]),as.data.frame(effect+residual))
if(orientation=="row") myY=cbind(NA,as.data.frame(effect+residual))

#Convert to category phenotype
if(category>0){
	myQuantile =(0:category)/category
	y.num= myY[,2]
	cutoff=quantile(y.num, myQuantile)
	y.cat= .bincode(y.num,cutoff,include.lowest = T)
	myY[,2]=y.cat
}

#print("Phenotype simulation accoplished")
return(list(Y=myY,u=effect,e=residual,QTN.position=QTN.position,effect=addeffect))
} #enf of function





