#Import library
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#Import data
setwd("~/Dropbox/gapit/GAPIT_Tutorial_Data")
myGD <- read.table("mdp_numeric.txt", head = TRUE)
myGM <- read.table("mdp_SNP_information.txt" , head = TRUE)
myG <- read.delim("mdp_genotype_test.hmp.txt", head = FALSE)

#Simulate phentype
set.seed(99164)
myPheno=GAPIT.Phenotype.Simulation(GD=myGD,h2=.5,NQTN=10,QTNDist="normal",effectunit=.299)
myY=myPheno$Y
colnames(myY)=c("taxa", "SimPheno")

#Override with modified functions
source('~/Dropbox/gapit/Functions/GAPIT.R', chdir = TRUE)
source('~/Dropbox/gapit/Functions/GAPIT.Main.R', chdir = TRUE)
source('~/Dropbox/gapit/Functions/GAPIT.Manhattan.R', chdir = TRUE)
source('~/Dropbox/gapit/Functions/GAPIT.MAF.R', chdir = TRUE)

setwd("~/Desktop/temp")
#GAPIT on G
myGAPIT <- GAPIT(
Y=myY,
G=myG,
group.from=1,
group.to=1,
group.by=10,
PCA.total=3,
memo="Test"
)

#GAPIT on GD
setwd("~/Desktop/temp")
myGAPIT=GAPIT(
plot.style="Beach",#options: Rainbow, Rushville, FarmCPU, Congress, Ocean, PLINK, Beach
Y=myY,
GD=myGD,
GM=myGM,
QTN.position= myPheno$QTN.position,
PCA.total=3,
group.from=1,
group.to=1,
group.by=10,
file.out=T,
#memo="Beach"
)

#BLINK
setwd("~/Desktop/temp")
myGD1=t(myGD[,-1])
write.table(myY,file="myData.txt",quote=F,row.name=F,col.name=T,sep="\t")
write.table(myGD1,file="myData.dat",quote=F,row.name=F,col.name=F,sep="\t")
write.table(myGM,file="myData.map",quote=F,row.name=F,col.name=T,sep="\t")
#run blink
system("~/Dropbox/zhanglab/BLINK/mac/blink --file myData --max_loop 10 --numeric --gwas")
#Extract p values
result<- read.table("SimPheno_GWAS_result.txt", head = TRUE)
myP=as.numeric(result[,5])
myGI.MP=cbind(myGM[,-1],myP)
GAPIT.Manhattan(GI.MP=myGI.MP,seqQTN=myPheno$QTN.position)
GAPIT.QQ(myP)

