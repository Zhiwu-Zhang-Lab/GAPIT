##############################################################################################
GAPIT <- function(Y,G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,turnOnEMMAxP3D=TRUE,
                groupFrom=1 ,groupTo=1000000,groupBy=10,CA="average", KT='Mean',Kinship.Method="Loiselle",
               	ngrid = 100, llin = -10, ulim = 10, esp = 1e-10,
                dataPath=NULL,numFiles=1,GFile=NULL, GFileExt=NULL,GDFile=NULL, GMFile=NULL, GDFileExt=NULL,GMFileExt=NULL,
                MAF.Filter.Rate=0.05,FDR.Rate = 0.05,FDR.Filter.Rate=1,model="Add",numPCs=0){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Outputv : See GAPIT.Main
#Authors: Zhiwu Zhang
# Last update: may 20, 2011 
GAPIT.Version="1.23"
##############################################################################################
#library('MASS')
#library(multtest)
if (ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")

for (trait in 2: ncol(Y))  {
myGAPIT <- GAPIT.Main(Y=Y[,c(1,trait)],G=G,GD=GD,GM=GM,KI=KI,Z=Z,CV=CV,turnOnEMMAxP3D=turnOnEMMAxP3D,Kinship.Method=Kinship.Method,
				              groupFrom=groupFrom,groupTo=groupTo,groupBy=groupBy,CA=CA,KT=KT,name.of.trait = colnames(Y)[trait],
                        dataPath=dataPath,numFiles=numFiles,GFile=GFile,GFileExt=GFileExt,GDFile=GDFile, GMFile=GMFile, GDFileExt=GDFileExt,GMFileExt=GMFileExt, 
                        MAF.Filter.Rate = MAF.Filter.Rate,FDR.Rate = FDR.Rate,FDR.Filter.Rate=FDR.Filter.Rate,model=model,,numPCs=numPCs,GAPIT.Version=GAPIT.Version)  
rm(myGAPIT)
gc()
}
print("GAPIT accomplished successfully!")
return()
}  #end of GAPIT function


##############################################################################################
GAPIT.Main <- function(Y,G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,turnOnEMMAxP3D=TRUE,
                groupFrom=1000000 ,groupTo=1,groupBy=10,CA="average", KT='Mean',Kinship.Method=NULL,
               	ngrid = 100, llin = -10, ulim = 10, esp = 1e-10,
                dataPath=NULL,numFiles=NULL,GFile=NULL, GFileExt=NULL,GDFile=NULL, GMFile=NULL, GDFileExt=NULL,GMFileExt=NULL,
                MAF.Filter.Rate=0.05,FDR.Rate = 0.05,FDR.Filter.Rate=1,model="Add",numPCs=0,  GAPIT.Version=GAPIT.Version,
                name.of.trait){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Output: GWAS table (text file), QQ plot (PDF), Manhattan plot (PDF), genomic prediction (text file), and 
#        genetic and residual variance components
#Authors: Zhiwu Zhang
# Last update: may 12, 2011 
##############################################################################################

seed=123
ratio=1

if(is.null(dataPath) & (is.null(G) & is.null(GD)) & (!is.null(GFile) | !is.null(GDFile))) stop("A path for genotype data should be provided!")
if(is.null(dataPath) & numFiles>1) stop(paste("A path for genotype data should be provided for the ",numFiles, " files specified by numFiles!"))

flag.single.vs.multiple=FALSE
if(!is.null(G)&numFiles>1) flag.single.vs.multiple=TRUE
if(!is.null(GD)&numFiles>1) flag.single.vs.multiple=TRUE
if(flag.single.vs.multiple) stop("Please input genotype data in either a singel file or use multiple file option!")

Timmer=GAPIT.Timmer(Infor="GAPIT (GWAS and GS in R)")
Memory=GAPIT.Memory(Infor="GAPIT (GWAS and GS in R)")

print(paste("Processing trait: ", name.of.trait)  )

#Multiple genotype files
if(is.null(numFiles))numFiles=1
if(numFiles==0)numFiles=1

hasGenotype=TRUE

print("Debug: reading G for first file")	
if(is.null(G) & !is.null(GFile) ) G <- read.table(paste(dataPath,GFile, "1.",GFileExt,sep=""), head = FALSE)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype loaded)")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype loaded")


if(!is.null(G)) genoFormat="Hapmap"
if(is.null(G)) genoFormat="EMMA"

#Rename GM as GI
if(!is.null(GM))GI=GM

if(is.null(GD) & !is.null(GDFile)  & !is.null(GMFile))
{
GI <- read.table(paste(dataPath,GMFile, "1.",GMFileExt,sep=""), head = TRUE)
GD <- read.table(paste(dataPath,GDFile, "1.",GDFileExt,sep=""), head = TRUE)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype data loaded)")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype data loaded")  
}

if(!is.null(GD) )
{
GT=as.matrix(GD[,1])  #get taxa
GD=GD[,-1] #remove taxa column
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GT created from GD)")
Memory=GAPIT.Memory(Memory=Memory,Infor="GT created from GD")  
}

if(!is.null(G) & !is.null(GD)) stop("Both hapmap and EMMA format exist, choose one only.")

if(!is.null(G))
{
#Convert HapMap to numerical
print(paste("Converting genotype...",sep=""))
hm=GAPIT.HapMap(G,model=model)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="HapMap")
Memory=GAPIT.Memory(Memory=Memory,Infor="HapMap")
print(paste("Converting genotype done.",sep=""))
rm(G)
gc()
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="G removed")
Memory=GAPIT.Memory(Memory=Memory,Infor="G removed")

GT=hm$GT
GD=hm$GD
GI=hm$GI
rm(hm)
gc()
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="hm removed")
Memory=GAPIT.Memory(Memory=Memory,Infor="hm removed")
}


#When GT and GD are missing, force to have fake ones (creating them from Y)
if(is.null(GD) & is.null(GT)) {
  hasGenotype=FALSE
	GT=data.frame(Y[,1])
	GD=matrix(1,nrow(Y),1)
}



#Handler of PCA and Kinship input on multiple geneotype files
if(is.null(KI) | numPCs>0)
GD.Sample=GAPIT.Sampler(GD=GD,dataPath=dataPath,numFiles=numFiles,GFile=GFile,GFileExt=GFileExt,seed=seed,ratio=ratio,
                        model=model, genoFormat=genoFormat, GDFile=GDFile, GDFileExt=GDFileExt)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Sampling genotype")
Memory=GAPIT.Memory(Memory=Memory,Infor="Sampling genotype")




#Examine status of datasets
if(is.null(Y)) stop ("Phenotype must exist.")
if(is.null(KI)&missing(GD)) stop ("Kinship is required. As genotype is not provided, kinship can not be created.")

#Create kinship from genotype if not provide
if(is.null(KI)&!is.null(GD)) 
{
print(dim(GD))


if(Kinship.Method=="EMMA")theKin= emma.kinship(snps=t(as.matrix(.5*GD.Sample)), method="additive", use="all")
if(Kinship.Method=="Loiselle")theKin= GAPIT.kinship.loiselle(snps=t(as.matrix(.5*GD.Sample)), method="additive", use="all")

KI=cbind(GT,as.data.frame(theKin))
rm(theKin)
gc()

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Estimating kinship")
Memory=GAPIT.Memory(Memory=Memory,Infor="Estimating kinship")

}


#Create Z as identity matrix from Y if it is not provided
if(is.null(Z)){
taxa=as.character(Y[,1])
Z=as.data.frame(diag(1,nrow(Y)))
Z=rbind(taxa,Z)
taxa=c('Taxa',as.character(taxa))
Z=cbind(taxa,Z)
}

#Add the part of non proportion in Z matrix
if(!is.null(Z)&nrow(Z)-1<nrow(Y)){
Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)
}

if(numPCs>0&!is.null(CV)){
CV=GAPIT.CVMergePC(CV,GAPIT.PCA(X = GD, taxa = GT, PC.number = numPCs))
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PCA")
Memory=GAPIT.Memory(Memory=Memory,Infor="PCA")
}


if(numPCs>0&is.null(CV)){
CV=GAPIT.PCA(X = GD.Sample, taxa = GT, PC.number = numPCs)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="EPCA")
Memory=GAPIT.Memory(Memory=Memory,Infor="PCA")
}




#Create CV with all 1's if it is not provided
if(is.null(CV)){
CV=Y[,1:2]
CV[,2]=1
colnames(CV)=c("taxa","overall")
}

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="ZMatrix")
#Memory=GAPIT.Memory(Memory=Memory,Infor="ZMatrix")

#Data quality control
qc <- GAPIT.QC(Y=Y,KI=KI, GT=GT,CV=CV,Z=Z)
GTindex=qc$GTindex
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="QC")
Memory=GAPIT.Memory(Memory=Memory,Infor="QC")

#Hanler of group boundry
if(groupFrom>groupTo) stop("groupTo should  be larger than groupFrom. Please correct them!")

if(!is.null(CV)& groupTo<ncol(CV)) {
#The minimum of group is number of columns in CV
  groupFrom=1 
  groupTo=1
  warning("The upper bound of groups (groupTo) is not sufficient. both boundries were set to a and GLM is performed!")
}

if(!is.null(CV)& groupFrom<1) {
  groupFrom=1 #minimum of group is number of columns in CV
  warning("The lower bound of groups should be 1 at least. Tt was set to 1!")
}

if(groupTo>nrow(qc$KI)) {
  groupTo=nrow(qc$KI) #maximum of group is number of rows in KI
	if(groupFrom>nrow(qc$KI)) groupFrom=nrow(qc$KI)
  warning("The upper bound of groups is too high. It was set to the size of kinship!")
}

#Optimization for group number, cluster algorithm and kinship type
GROUP=seq(groupTo,groupFrom,by=-groupBy)#The reverse order is to make sure to include full model 
if(missing("CA")) CA=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
if(missing("KT")) KT=c("Mean", "Max", "Min", "Median")
numSetting=length(GROUP)*length(CA)*length(KT)

optOnly=FALSE
if(numSetting>1 & hasGenotype) optOnly=TRUE

#Reform Y, GD and CV into EMMA format
ys=as.matrix(qc$Y[2])

X0=as.matrix(qc$CV[,-1])
xs.full=GD[qc$GTindex,]
xs.full=as.data.frame(xs.full)
xs.reduced=as.data.frame(xs.full[1,])

if( numSetting==1)
{ 
	xs=as.matrix(xs.full)
}else{
	xs=as.matrix(xs.reduced)
}
#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 4")
#Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 4")

#Initial
REMLs.best=10E20
ca.best=CA[1]
group.best=GROUP[1]
kt.best=KT[1]
count=0
compresion.KT<- list(NULL)
compresion.CA<- list(NULL)
compresion.GROUP<- list(NULL)
compresion.REML<- list(NULL)
compresion.VG<- list(NULL)
compresion.VE<- list(NULL)

#add intersection
overall <- rep(1, length(ys))
if(min(X0[,1])!=max(X0[,1])) X0 <- cbind(overall, X0) #do not add overall mean if X0 has it already at first column

#--------------------------------------------------------------------------------------------------------------------#



print("Compressing..." )
print(paste("The total combinations: ", numSetting,sep=""))
#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="DataProcessing")
#Memory=GAPIT.Memory(Memory=Memory,Infor="DataProcessing")

#Loop to optimize cluster algorithm, group number and kinship type 
for (ca in CA){
for (group in GROUP){
for (kt in KT){

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 1")
#Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 1")
count=count+1
          
if(group<ncol(X0)+1) group=1 # the emma function (emma.delta.REML.dLL.w.Z) does not allow K has dim less then CV

cp <- GAPIT.Compress(KI=qc$KI,CA=ca,KT=kt,GN=group,Timmer=Timmer,Memory=Memory)
Timmer=cp$Timmer
Memory=cp$Memory

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_cp")
Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2_cp")

bk <- GAPIT.Block(Z=qc$Z,GA=cp$GA,KG=cp$KG)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_bk")
Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 bk")

zc <- GAPIT.ZmatrixCompress(Z=qc$Z,GAU =bk$GA)

#write.table(zc$Z, "debug.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_zc")
Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 zc")

#Reform KW and Z into EMMA format
K=as.matrix(bk$KW)
z0=as.matrix(zc$Z[,-1])
Z=matrix(as.numeric(z0),nrow=nrow(z0),ncol=ncol(z0))

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 3")
#Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 3")



#Output kinship with phenotyped individuals for SAS
#taxa.Y=qc$Y[,1]
#taxa.K=qc$KI[,1]
#index=match(taxa.Y, taxa.K, nomatch = 0)
#KOnly=qc$KI[,-1]
#KS=as.matrix(KOnly[index,index])
#taxa.KS=taxa.K[index]
#order=order(taxa.KS)
#KSN=data.frame(cbind(as.character(taxa.KS[order]),KS))
#write.table(qc$Y, "Y.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(bk$KW, "KW.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(qc$CV, "CV.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(bk$GAU, "GAU.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(qc$GD, "GD.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(KSN, "KSN.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Prio PreP3D")
Memory=GAPIT.Memory(Memory=Memory,Infor="Prio PreP3D")

p3d <- GAPIT.EMMAxP3D(ys=ys,xs=xs,K = K ,Z=Z,X0=X0,GI=GI,turnOnEMMAxP3D=turnOnEMMAxP3D,Timmer=Timmer,Memory=Memory,
			 dataPath=dataPath,numFiles=numFiles,GFile=GFile,GFileExt=GFileExt,GDFile=GDFile, GMFile=GMFile, GDFileExt=GDFileExt,GMFileExt=GMFileExt,
       GTindex=qc$GTindex,genoFormat=genoFormat,optOnly=optOnly,model=model)
	
Timmer=p3d$Timmer
Memory=p3d$Memory

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Post PreP3D")
Memory=GAPIT.Memory(Memory=Memory,Infor="Post PreP3D")

#print("Cluster algorithm, kinship type, groups, VG, Ve and REML:")
print(paste(count, "of",numSetting,"--","Vg=",round(p3d$vgs,4), "VE=",round(p3d$ves,4),"-2LL=",round(p3d$REMLs,2), "  Clustering=",ca,"  Group number=", group ,"  Group kinship=",kt,sep = " "))

compresion.KT[[count]]=kt
compresion.CA[[count]]=ca
compresion.GROUP[[count]]=group
compresion.REML[[count]]=p3d$REMLs
compresion.VG[[count]]=p3d$vgs
compresion.VE[[count]]=p3d$ves


#Recording the best likelihood
if(p3d$REMLs!="NaN" || numSetting==1)
{
  if(p3d$REMLs<REMLs.best || numSetting==1)
  {
    REMLs.best=p3d$REMLs
    ca.best=ca
    group.best=group
    kt.best=kt
    cp.best=cp
    bk.best=bk
    zc.best=zc
    p3d.best=p3d
	
  }
}


}#end of for (ca in CA)

#Skip the rest group in case group 1 is finished
if(group==1) break #To skip the rest group iterations

}#end of for (group in GROUP)
}#end of for (kt in KT)
print(paste("Optimum: ",ca.best,kt.best,group.best,p3d.best$vgs, p3d.best$ves,p3d.best$REMLs ,sep = " "))
#--------------------------------------------------------------------------------------------------------------------#
#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Compression")
#Memory=GAPIT.Memory(Memory=Memory,Infor="Copmression")

if(numSetting==1) 
{
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS")
Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS")
}
#Collecting compression information
Compression= cbind(compresion.KT,compresion.CA,compresion.GROUP,compresion.REML,compresion.VG,compresion.VE)

#Perform GWAS with the optimum setting
#This section is omited if there is only one setting
print("Genomic screening..." )
if (numSetting>1) 
{

#Reform KW and Z into EMMA format
K=as.matrix(bk.best$KW)
z0=as.matrix(zc.best$Z[,-1])
Z=matrix(as.numeric(z0),nrow=nrow(z0),ncol=ncol(z0))

xs=as.matrix(xs.full)
p3d <- GAPIT.EMMAxP3D(ys=ys,xs=xs,K = K ,Z=Z,X0=X0,GI=GI,turnOnEMMAxP3D=turnOnEMMAxP3D,Timmer=Timmer,Memory=Memory,
	dataPath=dataPath,numFiles=numFiles,GFile=GFile,GFileExt=GFileExt,GDFile=GDFile, GMFile=GMFile, GDFileExt=GDFileExt,GMFileExt=GMFileExt,
  GTindex=qc$GTindex,genoFormat=genoFormat,optOnly=FALSE,model=model)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS")
Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS")

}

#Merge p value for multiple genotype files
if(numFiles>1)
{
theP <- read.table(paste("GAPIT.TMP.ps.","1",".txt",sep=""), head = FALSE)
theMAF <- read.table(paste("GAPIT.TMP.maf.","1",".txt",sep=""), head = FALSE)
thenobs <- read.table(paste("GAPIT.TMP.nobs.","1",".txt",sep=""), head = FALSE)
colnames(theP)="P"
colnames(theMAF )="MAF"
colnames(thenobs )="nobs"
ps=theP
MAF=theMAF
nobs=thenobs
for (file in 2:numFiles)
{
theGI <- read.table(paste("GAPIT.TMP.GI.",file,".txt",sep=""), head = TRUE)
theP <- read.table(paste("GAPIT.TMP.ps.",file,".txt",sep=""), head = FALSE)
theMAF <- read.table(paste("GAPIT.TMP.maf.",file,".txt",sep=""), head = FALSE)
thenobs <- read.table(paste("GAPIT.TMP.nobs.",file,".txt",sep=""), head = FALSE)
colnames(theP)="P"
colnames(theMAF )="MAF"
colnames(thenobs )="nobs"
ps=rbind(ps,theP)
MAF=rbind(MAF,theMAF)
nobs=rbind(nobs,thenobs)
GI=rbind(GI,data.frame(as.matrix(theGI)))


}
p3d$ps=ps
p3d$maf=MAF
p3d$nobs=nobs
}

#genomic prediction
if(length(bk.best$KW)>ncol(X0)) gs <- GAPIT.GS(KW=bk.best$KW,KO=bk.best$KO,KWO=bk.best$KWO,GAU=bk.best$GAU,UW=cbind(p3d.best$BLUP,p3d.best$PEV))
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GPS")
Memory=GAPIT.Memory(Memory=Memory,Infor="GPS")
 
#--------------------------------------------------------------------------------------------------------------------#
#Final report
print("Generating summary" )

#output BLUP and PEV
print("Genomic Breeding Values (GBV) ..." )
file=paste("GAPIT.", name.of.trait,".BLUP.csv" ,sep = "")
if(length(bk.best$KW)>ncol(X0)) write.table(gs$BLUP, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Make heatmap for distribution of BLUP and PEV
print("GBV and accuracy distribution..." )
if(length(bk.best$KW)>ncol(X0)) GAPIT.GS.Visualization(gsBLUP = gs$BLUP, BINS=20,name.of.trait = name.of.trait)

#Make a plot Summarzing the Compression Results, if more than one "compression level" has been assessed
print("Compression portfolios..." )
GAPIT.Compression.Visualization(Compression = Compression, name.of.trait = name.of.trait)

#Export GWAS results
if(hasGenotype ) 
{

print("Filtering SNPs with MAF..." )
PWI.Filtered <- GAPIT.Filter.SNPs.with.Low.MAFs(GI = GI, P.values = as.matrix(p3d$ps), maf = as.matrix(p3d$maf),
		     nobs = as.matrix(p3d$nobs), MAF.Filter.Rate = MAF.Filter.Rate)
  
#Run the BH multiple correction procedure of the results
#Create PWIP, which is a table of SNP Names, Chromosome, bp Position, Raw P-values, FDR Adjusted P-values
print("Calculating FDR..." )
PWIP <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered, FDR.Rate = FDR.Rate, FDR.Procedure = "BH")

  if(!is.null(PWI.Filtered)) 
  {
  #QQ plots
  print("QQ plot..." )
  GAPIT.QQ(P.values = PWIP$PWIP[,4], name.of.trait = name.of.trait,numSlots=10000)
  
  #Manhattan Plots
  print("Manhattan plot..." )
  GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, numSlots=10000, plot.type = "Chromosomewise")
  GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, numSlots=50000, plot.type = "Genomewise") 
  
  #Summary Table
  print("Association table..." )
  GAPIT.Table(final.table = PWIP$PWIP, name.of.trait = name.of.trait,FDR.Filter.Rate=FDR.Filter.Rate)
  }
}

#Log
log=GAPIT.Log(Y=Y,KI=KI,Z=Z,CV=CV,turnOnEMMAxP3D=turnOnEMMAxP3D,
				groupFrom = groupFrom ,groupTo =groupTo ,groupBy = groupBy ,CA = CA, KT= KT,
                      	ngrid = ngrid , llin = llin , ulim = ulim , esp = esp ,name.of.trait = name.of.trait)
#Memory usage
#GAPIT.Memory.Object(name.of.trait=name.of.trait)

#Timming
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Report")
Memory=GAPIT.Memory(Memory=Memory,Infor="Report")

file=paste("GAPIT.", name.of.trait,".Timming.csv" ,sep = "")
write.table(Timmer, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

file=paste("GAPIT.", name.of.trait,".Memory.Stage.csv" ,sep = "")
write.table(Memory, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

print(paste(name.of.trait, "has been analyzed successfully!") )
print(paste("The reulsts are saved in the directory of ", getwd()) )
print("---------------------------------------------------------------------------------------")

return (list(Timmer=Timmer,Compression=Compression))
}#The function GAPIT.Main ends here



##############################################################################################
GAPIT.EMMAxP3D <- function(ys,xs,K=NULL,Z=NULL,X0=NULL,GI=NULL,turnOnEMMAxP3D=TRUE,Timmer,Memory,ngrids=100,llim=-10,ulim=10,esp=1e-10,
		dataPath=NULL,numFiles=1,GFile=NULL,GFileExt=NULL,GTindex=NULL,GDFile=NULL, GMFile=NULL, GDFileExt=NULL,GMFileExt=NULL,genoFormat="Hapmap",
    optOnly=TRUE,model="Add"){
#Object: To esimate variance component by using EMMA algorithm and perform GWAS with P3D/EMMAx
#Output: ps, REMLs, stats, dfs, vgs, ves, BLUP,  BLUP_Plus_Mean, PEV
#Authors: Feng Tian, Alex Lipka and Zhiwu Zhang
# Last update: April 26, 2011 
# Library used: EMMA (Kang et al, Genetics, Vol. 178, 1709-1723, March 2008)
# Note: This function was modified from the function of emma.REML.t from the library
##############################################################################################

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="P3D Start")  
Memory=GAPIT.Memory(Memory=Memory,Infor="P3D Start")

#--------------------------------------------------------------------------------------------------------------------<
#Change data to matrix format if they are not
if (is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if (is.null(X0)) X0 <- matrix(1, ncol(ys), 1)

#handler of special Z and K
if(!is.null(Z)){ if(ncol(Z) == nrow(Z)) Z = NULL }  
if(!is.null(K)) {if(length(K)<2) K = NULL}

#Extract dimension information
g <- nrow(ys) #number of traits
n <- ncol(ys) #number of observation

q0 <- ncol(X0)#number of fixed effects
q1 <- q0 + 1  #Nuber of fixed effect including SNP  

nr=n
if(!is.null(K)) tv=ncol(K)

#decomposation without fixed effect
if(!is.null(K)) eig.L <- emma.eigen.L(Z, K) #this function handle both NULL Z and non-NULL Z matrix
if(!is.null(K)) eig.L$values[which(eig.L$values<0)]=0  #Negative eigen values
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.L")          
Memory=GAPIT.Memory(Memory=Memory,Infor="eig.L")



#decomposation with fixed effect (SNP not included)
X <-  X0 #covariate variables such as population structure
if (!is.null(Z) & !is.null(K)) eig.R <- emma.eigen.R.w.Z(Z, K, X) #This will be used to get REstricted ML (REML)
if (is.null(Z)  & !is.null(K)) eig.R <- emma.eigen.R.wo.Z(   K, X) #This will be used to get REstricted ML (REML)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.R")  
Memory=GAPIT.Memory(Memory=Memory,Infor="eig.R")




  
#-------------------------------------------------------------------------------------------------------------------->

#Loop on Traits
for (j in 1:g)
{

#--------------------------------------------------------------------------------------------------------------------<
if(!is.null(K)) REMLE <- emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R) 
if (!is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values + REMLE$delta)),rep(sqrt(1/REMLE$delta),nr - tv)),nr,((nr-tv)+length(eig.L$values)),byrow=TRUE)    
if ( is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(  sqrt(1/(eig.L$values + REMLE$delta)),nr,length(eig.L$values),byrow=TRUE)
x.prev <- vector(length = 0)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="REML")  
Memory=GAPIT.Memory(Memory=Memory,Infor="REML")
#-------------------------------------------------------------------------------------------------------------------->

if(optOnly) numFiles=1 #Skip sreening multiple genotype files

#Add loop for genotype data files
for (file in 1:numFiles)
{
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="New Genotype file")  
Memory=GAPIT.Memory(Memory=Memory,Infor="New Genotype file")

#update xs for each file
if(file>1)
{
rm(xs)
gc()

#Hapmap format
if(genoFormat=="Hapmap")
{
print("debug: reading G")
G <- read.table(paste(dataPath,GFile,file, ".",GFileExt,sep=""), head = FALSE)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype loaded")  
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype loaded")

#Convert HapMap to numerical
print(paste("Converting genotype...",sep=""))
hm=GAPIT.HapMap(G,model=model)
print(paste("Converting genotype done.",sep=""))
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="HapMap")  
Memory=GAPIT.Memory(Memory=Memory,Infor="HapMap")

G=NULL
gc()
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="G removed")  
Memory=GAPIT.Memory(Memory=Memory,Infor="G removed")

GT=hm$GT
print(paste("GT Generated.",sep=""))
Memory=GAPIT.Memory(Memory=Memory,Infor="GT Generated")

GI=hm$GI
print(paste("GI Generated.",sep=""))
Memory=GAPIT.Memory(Memory=Memory,Infor="GI Generated")

xs=hm$GD[GTindex,]
print(paste("xs Generated.",sep=""))
Memory=GAPIT.Memory(Memory=Memory,Infor="xs Generated")

rm(hm)
gc()
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="hm removed")  
Memory=GAPIT.Memory(Memory=Memory,Infor="hm removed")
} #nd of exising G

#EMMA format
if(genoFormat=="EMMA")
{
print("debug: reading GD and GM")
GI <- read.table(paste(dataPath,GMFile,file, ".",GMFileExt,sep=""), head = TRUE)
GD <- read.table(paste(dataPath,GDFile,file, ".",GDFileExt,sep=""), head = TRUE)
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype data loaded)")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype data loaded")  
GT=as.matrix(GD[,1])  #get taxa
GD=GD[,-1] #remove taxa column
xs=GD[GTindex,]
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GT created from GD)")
Memory=GAPIT.Memory(Memory=Memory,Infor="GT created from GD")  
}

} # end of multiple genotype files


if (is.null(dim(xs)) || nrow(xs) == 1)  xs <- matrix(xs, length(xs),1)
m <- ncol(xs) #number of SNPs
t <- nrow(xs) #number of individuals    

#allocate spaces for SNPs 
dfs <- matrix(nrow = m, ncol = g)
stats <- matrix(nrow = m, ncol = g)
ps <- matrix(nrow = m, ncol = g)
nobs <- matrix(nrow = m, ncol = g)
maf <- matrix(nrow = m, ncol = g)
#print(paste("Memory allocated.",sep=""))
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Memory allocation")  
Memory=GAPIT.Memory(Memory=Memory,Infor="Memory allocation")

if(optOnly)mloop=0
if(!optOnly)mloop=m

#Loop on SNPs
print(paste("Number of SNPs is ",mloop," in genotype file ",file, sep=""))
for (i in 0:mloop){
#print(i)
#--------------------------------------------------------------------------------------------------------------------<
    normalCase=TRUE  
    if ((i>0)&(floor(i/1000)==i/1000))  print(paste("Genotype file: ", file,", SNP: ",i," ",sep=""))

    # To extract current snp. It save computation for next one in case they are identical
    if(i==0&file==1){
      #For the model without fitting SNP
      vids <- !is.na(ys[j,]) #### Feng changed
      xv <- ys[j, vids]*0+1      #### Feng changed 
    }
    
    if(i>0){
      vids <- !is.na(xs[,i]) #### Feng changed
      xv <- xs[vids,i]      #### Feng changed 
      vids.TRUE=which(vids==TRUE)       
      vids.FALSE=which(vids==FALSE) 
      ns=length(xv)
      ss=sum(xv)
      maf[i]=min(.5*ss/ns,1-.5*ss/ns)
	nobs[i]=ns
    } 
  
    #Situation of no variation for SNP except the fisrt one(synthetic for EMMAx/P3D)
    if ((min(xv) ==max(xv) )&i>0) 
    {
      dfs[i, ] <- rep(NA, g)
      stats[i, ] <- rep(NA, g)
      ps[i, ] = rep(1, g)
      normalCase=FALSE      
    }

	 #Situation of the SNP is identical to previous    
    if (identical(x.prev, xv)) 
    { 
      dfs[i, ] <- dfs[i - 1, ]
      stats[i, ] <- stats[i - 1, ]
      ps[i, ] <- ps[i - 1, ]
      normalCase=FALSE       
    } 
#-------------------------------------------------------------------------------------------------------------------->

    
    #Normal case    
    if(normalCase)
    {  

#--------------------------------------------------------------------------------------------------------------------<
      #nv <- sum(vids)
      yv <- ys[j, vids] #### Feng changed
      nr <- sum(vids) #### Feng changed      
      if (!is.null(Z) & !is.null(K)) 
      {
        r<- ncol(Z) ####Feng, add a variable to indicate the number of random effect
        vran <- vids[1:r] ###Feng, add a variable to indicate random effects with nonmissing genotype
        tv <- sum(vran)  #### Feng changed
      }


      #Recalculate eig and REML if not using P3D
      if(turnOnEMMAxP3D==FALSE & !is.null(K)) 
      {
        if (!is.null(Z)) eig.R <- emma.eigen.R.w.Z(Z, K, X) #This will be used to get REstricted ML (REML)
        if (is.null(Z)) eig.R <- emma.eigen.R.wo.Z(   K, X) #This will be used to get REstricted ML (REML)
        if (!is.null(Z)) REMLE <- emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R) 
        if ( is.null(Z)) REMLE <- emma.REMLE(ys[j,], X, K, Z = NULL, ngrids, llim, ulim, esp, eig.R) 
        if (!is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values + REMLE$delta)),rep(sqrt(1/REMLE$delta),nr - tv)),nr,((nr-tv)+length(eig.L$values)),byrow=TRUE)    
        if ( is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(  sqrt(1/(eig.L$values + REMLE$delta)),nr,length(eig.L$values),byrow=TRUE)        
      }

  
   

#-------------------------------------------------------------------------------------------------------------------->

#--------------------------------------------------------------------------------------------------------------------<



      if(i>0) dfs[i, j] <- nr - q1
    	if(i>0) X <- cbind(X0[vids, , drop = FALSE], xs[vids,i])  
    	
      if(n==nr)
      {
        if(!is.null(K))
        {
            yt <- crossprod(U, yv)
            Xt <- crossprod(U, X)
        }else{
        yt=yv
        Xt=X
        }
        XX=crossprod(Xt, Xt)
  
  
  
        if(XX[1,1] == "NaN")
        {
          Xt[which(Xt=="NaN")]=0
          yt[which(yt=="NaN")]=0
          XX=crossprod(Xt, Xt)
        }
        XY=crossprod(Xt, yt)
      }
      
      #Missing SNP
      if(n>nr)
      {
       UU=crossprod(U,U)
       A11=UU[vids.TRUE,vids.TRUE]
       A12=UU[vids.TRUE,vids.FALSE]
       A21=UU[vids.FALSE,vids.TRUE]
       A22=UU[vids.FALSE,vids.FALSE]
       A22i=try(solve(A22) )
       if(inherits(A22i, "try-error")) A22i <- ginv(A22) 
       
       F11=A11-A12%*%A22i%*%A21
       XX=crossprod(X,F11)%*%X
       XY=crossprod(X,F11)%*%yv       
      }
      iXX <- try(solve(XX) )
      if(inherits(iXX, "try-error")) iXX <- ginv(crossprod(Xt, Xt)) 
      beta <- iXX %*% XY

#-------------------------------------------------------------------------------------------------------------------->

#--------------------------------------------------------------------------------------------------------------------<
      if(i==0 &file==1  & !is.null(K)) 
      {
        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="RducedModel")
	  Memory=GAPIT.Memory(Memory=Memory,Infor="ReducdModel")

        #Calculate the blup gammahat = vgKZprimeVinv*(Y-Xbetahat)
        vgs <- REMLE$vg
        ves <- REMLE$ve
        REMLs <- REMLE$REML  
        
        XtimesBetaHat <- X %*% beta 
        YminusXtimesBetaHat <- ys[j,]- XtimesBetaHat
        vgK <- REMLE$vg*K
        Dt <- crossprod(U, YminusXtimesBetaHat)        
        if (!is.null(Z))  Zt <- crossprod(U, Z)   
        if (is.null(Z)) Zt <- t(U) 

        if(XX[1,1] == "NaN")
        {
        Dt[which(Dt=="NaN")]=0
        Zt[which(Zt=="NaN")]=0
        }
         
        BLUP <- K %*% crossprod(Zt, Dt) #Using K instead of vgK because using H=V/Vg

        grand.mean.vector <- rep(beta[1], length(BLUP))
        BLUP_Plus_Mean <- grand.mean.vector + BLUP 
    	  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="BLUP")          
        Memory=GAPIT.Memory(Memory=Memory,Infor="BLUP")

        #PEV
        C11=try(vgs*solve(crossprod(Xt,Xt)))
        if(inherits(C11, "try-error")) C11=vgs*ginv(crossprod(Xt,Xt))
        
        C21=-K%*%crossprod(Zt,Xt)%*%C11
        Kinv=try(solve(K)    )
        if(inherits(Kinv, "try-error")) Kinv=ginv(K)

        if(!is.null(Z)) term.0=crossprod(Z,Z)/ves
        if(is.null(Z)) term.0=diag(1/ves,nrow(K))
        
        term.1=try(solve(term.0+Kinv/vgs )  )
        if(inherits(term.1, "try-error")) term.1=ginv(term.0+Kinv/vgs )
        
        term.2=C21%*%crossprod(Xt,Zt)%*%K
        C22=(term.1-term.2 )
        PEV=as.matrix(diag(C22))
        BLUE=X%*%beta

    	  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PEV") 
        Memory=GAPIT.Memory(Memory=Memory,Infor="PEV")        
          
      }#end of if(i==0&file==1   & !is.null(K)) 
#-------------------------------------------------------------------------------------------------------------------->

#--------------------------------------------------------------------------------------------------------------------<
      if(i==0 &file==1  & is.null(K)) 
      {
        YY=crossprod(yt, yt)  
        ves=(YY-crossprod(beta,XY))/(n-q0)          
        r=yt-X%*%iXX%*%XY
        REMLs=-.5*(n-q0)*log(det(ves))                  -.5*n                   -.5*(n-q0)*log(2*pi)        
#       REMLs=-.5*n*log(det(ves)) -.5*log(det(iXX)/ves) -.5*crossprod(r,r)/ves  -.5*(n-q0)*log(2*pi)  
        vgs = 0
        BLUP = 0 
        BLUP_Plus_Mean = NaN
        PEV = ves
       	BLUE=X%*%beta
      }
      
      #calculate t statistics and probabilty 
      if(i > 0)	
      { 
        if(!is.null(K)) stats[i, j] <- beta[q1]/sqrt(iXX[q1, q1] *REMLE$vg)
        if(is.null(K)) stats[i, j] <- beta[q1]/sqrt(iXX[q1, q1] *ves)
        ps[i, ] <- 2 * pt(abs(stats[i, ]), dfs[i, ],lower.tail = FALSE)
      }
#-------------------------------------------------------------------------------------------------------------------->

    } # End of if(normalCase)
    x.prev=xv #update SNP

} # End of loop on SNPs  

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Screening SNPs")  
Memory=GAPIT.Memory(Memory=Memory,Infor="Screening SNPs")

#output p value for the genotype file
if(numFiles>1) 
{
  write.table(GI, paste("GAPIT.TMP.GI.",file,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(ps, paste("GAPIT.TMP.ps.",file,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  write.table(maf, paste("GAPIT.TMP.maf.",file,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  write.table(nobs, paste("GAPIT.TMP.nobs.",file,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  #rm(dfs,stats,ps,nobs,maf,GI)   #This cause problem on return
  #gc()
}
    
} # Ebd of loop on file 
} # End of loop on traits

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS done for this Trait")
Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS done for this Trait")


return(list(ps = ps, REMLs = -2*REMLs, stats = stats, dfs = dfs,maf=maf,nobs = nobs,Timmer=Timmer,Memory=Memory,
        vgs = vgs, ves = ves, BLUP = BLUP, BLUP_Plus_Mean = BLUP_Plus_Mean,
        PEV = PEV, BLUE=BLUE))

#print("GAPIT.EMMAxP3D accomplished successfully!")
}#end of GAPIT.EMMAxP3D function

##############################################################################################
GAPIT.ZmatrixFormation <- function(Z,Y){
#Object: To expande the proportion Z to final Z
#Output: Z
#Authors: Zhiwu Zhang 
# Last update: April 22, 2011 
##############################################################################################

#split individuals in Y to the ones that are given Z and the one not
taxa.Z=as.matrix(Z[-1,1])
taxa.Y=as.matrix(Y[,1])
taxa.diff=setdiff(taxa.Y,taxa.Z)
taxa.I=as.matrix(taxa.Y[match(taxa.diff,taxa.Y,nomatch = 0)])
taxa.Z.col=as.matrix(Z[1,-1])

#Create final Z with zero block and identity block
Z0=matrix(data=0,nrow=nrow(taxa.Z),ncol=nrow(taxa.I))
Z1=diag(1,nrow(taxa.I))
ZC=as.matrix(rbind(Z0,Z1))

#To label rows and columns
label.row=rbind(as.matrix(Z[,1]),taxa.I)
label.col=t(taxa.I)

#update the zero block by the given Z matrix
position=t(as.matrix(match(taxa.Z.col,taxa.I,nomatch = 0)))
ZC[1:nrow(taxa.Z),position]=as.matrix(Z[-1,-1])

#habdler of parents do not have phenotype (colums of Z are not in taxa.I)
# To do list

#To form final Z matrix
dataPart=rbind(label.col,ZC)
Z=data.frame(cbind(label.row,dataPart))

#print("GAPIT.ZmatrixFormation accomplished successfully!")
return(Z)
}#The function GAPIT.ZmatrixFormation ends here

##############################################################################################
GAPIT.QC <- function(Y,KI,GT,CV,Z){
#Object: to do data quality control
#Output: Y, KI, GD, CV, Z, flag
#Authors: Zhiwu Zhang and Alex Lipka 
# Last update: April 14, 2011 
##############################################################################################

#Remove missing phenotype
Y=Y[which(Y[,2]!="NaN"),]

# Remove duplicates: 
# GT row wise, Z column wise, and KI both direction.
if(exists("GT"))
{ 
taxa.kept=unique(GT[,1])
}

taxa.all=KI[,1]
taxa.uniqe=unique(taxa.all)
position=match(taxa.uniqe, taxa.all,nomatch = 0)
position.addition=cbind(1,t(1+position))
KI=KI[position,position.addition]

if(exists("Z"))
{
taxa.all=as.matrix(Z[1,])
taxa.uniqe=intersect(taxa.all,taxa.all)
position=match(taxa.uniqe, taxa.all,nomatch = 0)
Z=Z[,position]
}

#Remove the columns of Z if they are not in KI/GT. KI/GT are allowed to have individuals not in Z
taxa.all=KI[,1]
taxa.kinship=unique(taxa.all)
taxa.Z=as.matrix(Z[1,])
#taxa.Z=colnames(Z) #This does not work for names starting with numerical or "-"
taxa.Z_K_common=intersect(taxa.kinship,taxa.Z)
Z <-cbind(Z[,1], Z[,match(taxa.Z_K_common, taxa.Z, nomatch = 0)])

#Remove the rows of Z if all the ellements sum to 0
Z1=Z[-1,-1]
Z2=data.frame(Z1)
Z3=as.matrix(Z2)
Z4=as.numeric(Z3) #one dimemtion
Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
RS=rowSums(Z5)>0
#The above process could be simplified!
Z <- Z[c(TRUE,RS),]

#make individuals the same in Z, Y, GT and CV
# get intersect of all the data
taxa=intersect(Z[-1,1],Y[,1])
if(exists("GT"))taxa=intersect(taxa,taxa.kept)
if(exists("CV"))taxa=intersect(taxa,CV[,1])

#keep the common ones
t=c(TRUE, Z[-1,1]%in%taxa)
Z <- Z[t,]

#Remove the columns of Z if all the ellements sum to 0
Z1=Z[-1,-1]
Z2=data.frame(Z1)
Z3=as.matrix(Z2)
Z4=as.numeric(Z3) #one dimemtion
Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
CS=colSums(Z5)>0
#The above process could be simplified!
Z <- Z[,c(TRUE,CS)]


Y <- Y[Y[,1]%in%taxa,]
if(exists("GT")) taxa.kept=data.frame(taxa.kept[taxa.kept%in%taxa])
if(exists("CV")) CV=CV[CV[,1]%in%taxa,]

#get position of taxa.kept in GT
position=match(taxa.kept[,1], GT[,1],nomatch = 0)

#To sort Y, GT, CV and Z
Y=Y[order(Y[,1]),]
CV=CV[order(CV[,1]),]
order.taxa.kept=order(taxa.kept[,1])
Z=Z[c(1,1+order(Z[-1,1])),]
GTindex=position[order.taxa.kept]
flag=nrow(Y)==nrow(Z)-1&nrow(Y)==nrow(GT)&nrow(Y)==nrow(CV)

#print("GAPIT.QC accomplished successfully!")
return(list(Y = Y, KI = KI, GT = GT, CV = CV, Z = Z, GTindex=GTindex, flag=flag))
}#The function GAPIT.QC ends here

##############################################################################################
GAPIT.Compress <- function(KI,CA = "average",KT = "Mean",GN=nrow(KI),Timmer,Memory){
#Object: To cluster individuals into groups based on kinship
#Output: GA, KG
#Authors: Alex Lipka and Zhiwu Zhang 
# Last update: April 14, 2011 
##############################################################################################

#For debug 
#GN=7

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP start") 
Memory=GAPIT.Memory(Memory=Memory,Infor="cp start")

# Extract the line names
line.names <- KI[,1]

# Remove the first column of the kinship matrix, which is the line names
KI <- KI[ ,-1]

# Convert kinship to distance
distance.matrix <- 2 - KI 
distance.matrix.as.dist <- as.dist(distance.matrix)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP distance") 
Memory=GAPIT.Memory(Memory=Memory,Infor="cp distance")

# hclust() will perform the hiearchical cluster analysis
cluster.distance.matrix <- hclust(distance.matrix.as.dist, method = CA)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP cluster") 
Memory=GAPIT.Memory(Memory=Memory,Infor="cp cluster")

# Cutree will assign lines into k clusters
group.membership <- cutree(cluster.distance.matrix, k = GN) 

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP cutree") 
#Memory=GAPIT.Memory(Memory=Memory,Infor="cp cutree")

#calculate group kinship
if(KT == "Mean"){
#This matrix operation is much faster than tapply function for  "Mean"
x=as.factor(group.membership)
#b = model.matrix(~x-1) 
n=max(as.numeric(as.vector(x)))
b=diag(n)[x,]

KG=t(b)%*%as.matrix(KI)%*%b
CT=t(b)%*%(0*as.matrix(KI)+1)%*%b
KG=as.matrix(KG/CT)
rownames(KG)=c(1:nrow(KG))
colnames(KG)=c(1:ncol(KG))
#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP calculation original")
#Memory=GAPIT.Memory(Memory=Memory,Infor="cp calculation original")
 

}

gm=as.factor(group.membership)
kv=as.numeric(as.matrix(KI))
kvr=rep(gm,ncol(KI))
kvc=as.numeric(t(matrix(kvr,nrow(KI),ncol(KI))))
kInCol=t(rbind(kv,kvr,kvc))




if(KT != "Mean"){
#if (KT == "Mean")
#    KG<- tapply(kInCol[,1], list(kInCol[,2], kInCol[,3]), mean)

if (KT == "Max")    
    KG <- tapply(kInCol[,1], list(kInCol[,2], kInCol[,3]), max)
if (KT == "Min")   
    KG <- tapply(kInCol[,1], list(kInCol[,2], kInCol[,3]), min)    
if (KT == "Median")  
    KG <- tapply(kInCol[,1], list(kInCol[,2], kInCol[,3]), median)  
} #this is end of brancing "Mean" and the rest
    
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP calculation") 
Memory=GAPIT.Memory(Memory=Memory,Infor="cp calculation")

# add line names 
GA <- data.frame(group.membership)
GA <- data.frame(cbind(as.character(line.names),as.numeric(group.membership) ))

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP Final") 
#Memory=GAPIT.Memory(Memory=Memory,Infor="CP Final")

#print("GAPIT.Compress accomplished successfully!")
return(list(GA=GA, KG=KG,Timmer=Timmer,Memory=Memory))
}#The function GAPIT.Compress ends here

##############################################################################################
GAPIT.Block <- function(Z,GA,KG){
#Object: To split a group kinship into two blocks containing individuals with and without phenotype
#Output: GAU,KW,KO,KWO
#Authors: Zhiwu Zhang and Alex Lipka 
# Last update: April 14, 2011 
##############################################################################################

# To separate group kiship into two blocks: with and without phenotype.
# A group goes to with phenotype as loog as it has one phenotyped individual.

#find position in group assignment (GA) for the individual associate with phenotype (specified by Z)
#taxa=unique(intersect(as.matrix(Z[1,-1]),GA[,1]))

taxa.Z=as.matrix(Z[1,-1])
taxa.GA=as.matrix(GA[,1])
position=taxa.GA%in%taxa.Z

#Initial block as 2
GAU=cbind(GA,2)

#Assign block as 1 if the individual has phenotype
GAU[position,3]=1

#Modify the non-phenotyped individuals if they in a group with phenotyped individuals
#To find the groups with phenotyped individuals
#update block assignment for all these groups
#get list of group that should be block 1

grp.12=as.matrix(unique(GAU[,2]))
grp.1=as.matrix(unique(GAU[which(GAU[,3]==1),2]))
grp.2= as.matrix(setdiff(grp.12,grp.1))
numWithout=length(grp.2)

order.1=1:length(grp.1)
order.2=1:length(grp.2)
if(numWithout >0) grpblock=as.matrix(rbind(cbind(grp.1,1,order.1), cbind(grp.2,2,order.2)))
if(numWithout==0) grpblock=as.matrix(      cbind(grp.1,1,order.1),                       )

order.block=order(as.matrix(GAU[,3]))
colnames(grpblock)=c("grp","block","ID")

GAU0 <- merge(GAU[order.block,-3], grpblock, by.x = "X2", by.y = "grp")
GAU=GAU0[,c(2,1,3,4)]

KW=KG[grp.1,grp.1]
KO=KG[grp.2,grp.2]
KWO=KG[grp.1,grp.2]

#write.table(GAU, "GAU.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

#print("GAPIT.Block accomplished successfully!")
return(list(GAU=GAU,KW=KW,KO=KO,KWO=KWO))
}#The function GAPIT.Block ends here


##############################################################################################
GAPIT.ZmatrixCompress <- function(Z,GAU){
#Object: To assign the fraction of a individual belonging to a group
#Output: Z
#Authors: Zhiwu Zhang
# Last update: April 14, 2011 
##############################################################################################

#sort Z column wise
order.Z=order(as.matrix(Z[1,-1]))
Z1=Z[-1,-1]
Z1 <- Z1[,order.Z]

#Extraction of GAU coresponding to Z, sort GAU rowwise to mach columns of Z, and make design matrix
effect.Z=as.matrix(Z[1,-1])
effect.GAU=as.matrix(GAU[,1])
GAU0=GAU[effect.GAU%in%effect.Z,]
order.GAU=order(GAU0[,1])
GAU1 <- GAU0[order.GAU,]
id.1=GAU1[which(GAU1[,3]==1),4]
n=max(as.numeric(as.vector(id.1)))
x=as.numeric(as.matrix(GAU1[,4]))
DS=diag(n)[x,]

#write.table(b, "debug.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
#write.table(GAU1, "debug2.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

#Z matrix from individual to group
Z1.numeric <- as.numeric(as.matrix(Z1))
Z1.matrix <- matrix(Z1.numeric, nrow = nrow(Z1), ncol = ncol(Z1)) 
Z2=Z1.matrix%*%DS

Z3=data.frame(cbind(as.character(Z[-1,1]),Z2))
Z=Z3[order(Z3[,1]),]


#write.table(DS, "debug.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

#print("GAPIT.ZmatrixCompress accomplished successfully!")
return(list(Z=Z))
}#The function GAPIT.ZmatrixCompress ends here

##############################################################################################
GAPIT.GS <- function(KW,KO,KWO,GAU,UW){
#Object: to derive BLUP for the individuals without phenotype
#Output: BLUP
#Authors: Zhiwu Zhang 
# Last update: April 17, 2011 
##############################################################################################

UO=try(t(KWO)%*%solve(KW)%*%UW)
if(inherits(UO, "try-error")) UO=t(KWO)%*%ginv(KW)%*%UW

n=ncol(UW) #get number of columns, add additional for individual name

#Assign BLUP of group to its individuals
BLUP=data.frame(as.matrix(GAU[,1:4]))

BLUP.W=BLUP[which(GAU[,3]==1),]
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
#write.table(index.W, "debug.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

#print("GAPIT.GS accomplished successfully!")
return(list(BLUP=all))
}#The function GAPIT.GS ends here

##############################################################################################
GAPIT.Numericalization.Add <- function(x){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################
n=length(x)
lev=levels(as.factor(x))
lev=setdiff(lev,"NN")
len=length(lev)
#print(lev)

#x1=ifelse(x=="NN",NA,ifelse(x==lev[1],0,ifelse(x==lev[2],1,2)))

#1 or less status
if(len<=1)x=0

#2 status
if(len==2)x=ifelse(x=="NN",1,ifelse(x==lev[1],0,2))

#3 status
if(len==3)x=ifelse(x=="NN",1,ifelse(x==lev[1],0,ifelse(x==lev[2],1,2)))

# more than 3 status
if(len> 3){
x=0
warning("Existing none biallelic SNP!")
}

return(matrix(x,n,1))
}#end of GAPIT.Numericalization.Add function

##############################################################################################
GAPIT.Numericalization.Dom <- function(x){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################

lev=levels(as.factor(x))
lev=setdiff(lev,"NN")
len=length(lev)
#print(lev)

#x1=ifelse(x=="NN",NA,ifelse(x==lev[1],0,ifelse(x==lev[2],1,2)))

#1 or less status
if(len<=1)x1=0

# more than 3 status
if(len> 3){
x1=0
warning("Existing none biallelic SNP!")
}

#2 status
if(len==2|len==3)x1=ifelse(x=="NN",1,ifelse(substr(x,1,1)==substr(x,2,2),0,2))

return(matrix(x1,length(x),1))
}#end of GAPIT.Numericalization.Dom function

##############################################################################################
GAPIT.HapMap <- function(G,model="Add"){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################
gc()
GAPIT.Memory.Object(name.of.trait="HapMap.Start")

#GT=data.frame(G[1,-(1:11)])
GT= t(G[1,-(1:11)])
GI= G[-1,c(1,3,4)]
if(model=="Add")GD= apply(G[-1,-(1:11)],1,GAPIT.Numericalization.Add)
if(model!="Add")GD= apply(G[-1,-(1:11)],1,GAPIT.Numericalization.Dom)

colnames(GT)="taxa"
colnames(GI)=c("SNP","Chromosome","Position")

GAPIT.Memory.Object(name.of.trait="HapMap.Finished")
return(list(GT=GT,GD=GD,GI=GI))
}#end of GAPIT.HapMap function

##############################################################################################
GAPIT.CVMergePC<- function(X,Y){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################

#Z=X+Y

Z <- merge(X, Y, by.x = colnames(X)[1], by.y = colnames(Y)[1])

return(Z)
}#end of GAPIT.CVMergePCfunction



##############################################################################################
GAPIT.Memory <- function(Memory =NULL,Infor){
#Object: To report memory usage
#Output: Memory 
#Authors: Zhiwu Zhang
# Last update: June 6, 2011 
##############################################################################################
gc()
size <- memory.size()
#print(paste("Memory usage: ",size," for", Infor))
if(is.null(Memory)) {
Increased=0
Memory =cbind(Infor,size ,Increased)
}else{
Increased=0
Memory.current=cbind(Infor,size ,Increased)
Memory=rbind(Memory,Memory.current)
Memory[nrow(Memory),3]=as.numeric(as.matrix(Memory[nrow(Memory),2]))-as.numeric(as.matrix(Memory[nrow(Memory)-1,2]))
}

return (Memory)
}#end of GAPIT.Memory function

##############################################################################################
GAPIT.Memory.Object <- function(name.of.trait="Trait"){
# Object: To report memoery usage
# Authors: Heuristic Andrew
# http://heuristically.wordpress.com/2010/01/04/r-memory-usage-statistics-variable/
# Modified by Zhiwu Zhang
# Last update: may 29, 2011 
##############################################################################################
  
# print aggregate memory usage statistics 
print(paste('R is using', memory.size(), 'MB out of limit', memory.limit(), 'MB')) 
  
# create function to return matrix of memory consumption 
object.sizes <- function() 
{ 
    return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
        object.size(get(object.name)))))) 
} 

# export file in table format 
memory=object.sizes() 
file=paste("GAPIT.", name.of.trait,".Memory.Object.csv" ,sep = "")
write.table(memory, file, quote = FALSE, sep = ",", row.names = TRUE,col.names = TRUE)


# export file in PDF format 
pdf(paste("GAPIT.", name.of.trait,".Memory.Object.pdf" ,sep = ""))
# draw bar plot 
barplot(object.sizes(), 
    main="Memory usage by object", ylab="Bytes", xlab="Variable name", 
    col=heat.colors(length(object.sizes()))) 
# draw dot chart 
dotchart(object.sizes(), main="Memory usage by object", xlab="Bytes") 
# draw pie chart 
pie(object.sizes(), main="Memory usage by object")
dev.off()  
}  



##############################################################################################
GAPIT.Timmer <- function(Timmer=NULL,Infor){
#Object: To report current time
#Output: Timmer
#Authors: Zhiwu Zhang
# Last update: may 8, 2011 
##############################################################################################

Time<- Sys.time()
if(is.null(Timmer)) {
Elapsed=0
Timmer=cbind(Infor,Time,Elapsed)
}else{
Elapsed=0
Timmer.current=cbind(Infor,Time,Elapsed)
Timmer=rbind(Timmer,Timmer.current)
Timmer[nrow(Timmer),3]=as.numeric(as.matrix(Timmer[nrow(Timmer),2]))-as.numeric(as.matrix(Timmer[nrow(Timmer)-1,2]))
}

#print(paste('Time used: ', Timmer[nrow(Timmer),3], ' seconds for ',Infor,sep="" )) 
return (Timmer)
}#end of GAPIT.EMMAxP3D function


##############################################################################################
GAPIT.Log <- function(Y=Y,KI=KI,Z=Z,CV=CV,turnOnEMMAxP3D=turnOnEMMAxP3D,
				groupFrom = groupFrom ,groupTo =groupTo ,groupBy = groupBy ,CA = CA, KT= KT,
                      	ngrid = ngrid , llin = llin , ulim = ulim , esp = esp ,name.of.trait = name.of.trait){
#Object: To report model factors
#Output: Text file (GAPIT.Log.txt)
#Authors: Zhiwu Zhang
# Last update: may 16, 2011 
##############################################################################################

#Creat storage
facto <- list(NULL)
value <- list(NULL)

#collecting model factors

facto[[1]]="Trait"
value[[1]]=paste(dim(Y))

facto[[2]]="groupBy "
value[[2]]=groupBy 

facto[[3]]="Trait name "
value[[3]]=name.of.trait

facto[[4]]="Kinship"
value[[4]]=dim(KI)

facto[[5]]="Z Matrix"
value[[5]]=dim(Z)

facto[[6]]="Covariate"
value[[6]]=dim(CV)

facto[[7]]="EMMAxP3D"
value[[7]]=turnOnEMMAxP3D

facto[[8]]="Clustering algorithms"
value[[8]]=CA

facto[[9]]="Group kinship"
value[[9]]=KT

facto[[10]]="groupFrom "
value[[10]]=groupFrom 

facto[[11]]="groupTo "
value[[11]]=groupTo 



theLog=as.matrix(cbind(facto,value))
#theLog=as.character(as.matrix(cbind(facto,value)))
colnames(theLog)=c("Model", "Value")
file=paste("GAPIT.", name.of.trait,".Log.csv" ,sep = "")
write.table(theLog, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

return (theLog)
}


##############################################################################################
GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure <- function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
#Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
#Output: PWIP, number.of.significant.SNPs
#Authors: Alex Lipka and Zhiwu Zhang 
# Last update: May 5, 2011 
##############################################################################################


    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest

if(is.null(PWI))
{
PWIP=NULL
number.of.significant.SNPs = 0
}

if(!is.null(PWI))
{  
 
    #library(multtest)
    
    if(dim(PWI)[1] == 1){
     PWIP <- cbind(PWI, PWI[4])
     colnames(PWIP)[5] <- "FDR_Adjusted_P-values"
    }
   
    if(dim(PWI)[1] > 1){ 
    #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
    #Right column: Simes corrected p-value
    res <- mt.rawp2adjp(PWI[,4], FDR.Procedure)

    #This command should order the p-values in the order of the SNPs in the data set
  adjp <- res$adjp[order(res$index), ]

  #round(adjp[1:7,],4)
    #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
  temp <- mt.reject(adjp[,2], FDR.Rate)

    #Lists all number of SNPs that were rejected by the BY procedure
  #temp$r

    #Attach the FDR adjusted p-values to AS_Results

  PWIP <- cbind(PWI, adjp[,2])

    #Sort these data by lowest to highest FDR adjusted p-value
  PWIP <- PWIP[order(PWIP[,4]),]
  
  colnames(PWIP)[7] <- "FDR_Adjusted_P-values"
  number.of.significant.SNPs = temp$r
  }
  #print("GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure accomplished successfully!")
}  
  return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))

}



##############################################################################################
GAPIT.Pruning <- function(values,numSlots=20){
#Object: To get index of subset that evenly distribute
#Output: Index
#Authors: Zhiwu Zhang
# Last update: May 28, 2011 
##############################################################################################

#values= log.P.values

values=sqrt(values)  #This shift the weight a little bit to the low building.
theMin=min(values)
theMax=max(values)
range=theMax-theMin
interval=range/numSlots

ladder=round(values/interval)
ladder2=c(ladder[-1],0)
keep=ladder-ladder2
index=which(keep>0)

return(index)
}#end of GAPIT.Pruning 



##############################################################################################
GAPIT.QQ <- function(P.values, plot.type = "log_P_values", name.of.trait = "Trait",numSlots=1000){
#Object: Make a QQ-Plot of the P-values
#Options for plot.type = "log_P_values" and "P_values" 
#Output: A pdf of the QQ-plot
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: May 9, 2011 
##############################################################################################


# Sort the data by the raw P-values
P.values <- P.values[order(P.values)]
  
#Set up the p-value quantiles
p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)


if(plot.type == "log_P_values")
{
    log.P.values <- -log10(P.values)
    log.Quantiles <- -log10(p_value_quantiles)
	
    index=GAPIT.Pruning(log.P.values,numSlots=numSlots)
    log.P.values=log.P.values[index ]
    log.Quantiles=log.Quantiles[index]

    pdf(paste("GAPIT.", name.of.trait,".QQ-Plot.pdf" ,sep = ""))
    par(mar = c(5,5,5,5))
    qqplot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), 
           cex.axis=1.5, cex.lab=2, lty = 1, lwd = 1, col = "Blue" ,xlab =expression(Expected~~-log[10](italic(p))),
           ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=" "))
    abline(a = 0, b = 1, col = "red")
    dev.off()   
}


if(plot.type == "P_values")
{
  pdf(paste("QQ-Plot_", name.of.trait,".pdf" ,sep = ""))
  par(mar = c(5,5,5,5))
  qqplot(p_value_quantiles, P.values, xlim = c(0,1), 
         ylim = c(0,1), type = "l" , xlab = "Uniform[0,1] Theoretical Quantiles", 
         lty = 1, lwd = 1, ylab = "Quantiles of P-values from GWAS", col = "Blue",
         main = paste(name.of.trait,sep=" "))
  abline(a = 0, b = 1, col = "red")
  dev.off()   
}


  #print("GAPIT.QQ  accomplished successfully!")
  

}


##############################################################################################
GAPIT.Manhattan <- function(GI.MP = NULL, name.of.trait = "Trait", 
                   plot.type = "Genomewise",numSlots=1000){
#Object: Make a Manhattan Plot
#Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome" 
#Output: A pdf of the Manhattan Plot
#Authors: Alex Lipka, Zhiwu Zhang, and Meng Li 
# Last update: May 10, 2011 
##############################################################################################

#do nothing if null input
if(is.null(GI.MP)) return

GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))

#Remove all SNPs that do not have a choromosome and bp position
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]

#Remove all SNPs that have P values above 0
GI.MP <- GI.MP[GI.MP[,3]>0,]

#Replace P the -log10 of the P-values
GI.MP[,3] <-  -log10(GI.MP[,3])
y.lim <- ceiling(max(GI.MP[,3]))
chm.to.analyze <- unique(GI.MP[,1]) 
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)

#Chromosomewise plot
if(plot.type == "Chromosomewise")
{
  pdf(paste("GAPIT.", name.of.trait,".Manhattan-Plot.Chromosomewise.pdf" ,sep = ""), width = 10)
  par(mar = c(5,5,4,3), lab = c(8,5,7))
  for(i in 1:numCHR)
  {
    subset=GI.MP[GI.MP[,1]==chm.to.analyze[i],]
    x0 <- as.numeric(subset[,2])/10^(8)
    y0 <- as.numeric(subset[,3])
	order=order(y0,decreasing = TRUE)
    index0=GAPIT.Pruning(y0[order],numSlots=20000)
	index=order(index0)
	x=x0[index]
	y=y0[index]
    #color.vector <- subset(temp.par.data[,7], temp.par.data[,4] == i)
    plot(y~x,type="p", ylim=c(0,y.lim), xlim = c(min(x), max(x)), col = "navy", xlab = expression(Base~Pairs~(x10^-8)), ylab = "-Log Base 10 p-value", main = paste("Chm",chm.to.analyze[i],sep=" "))
  }
  dev.off()
} #Chromosomewise plot


#Genomewise plot
if(plot.type == "Genomewise")
{
  
  GI.MP <- GI.MP[order(GI.MP[,2]),]
  GI.MP <- GI.MP[order(GI.MP[,1]),]
  color.vector <- rep(c("orangered","navyblue"),numCHR)
  ticks=NULL
  lastbase=0
  
  #change base position to accumulatives
  for (i in chm.to.analyze)
  {
    index=(GI.MP[,1]==i)
    ticks <- c(ticks, lastbase+mean(GI.MP[index,2])) 
    GI.MP[index,2]=GI.MP[index,2]+lastbase
    lastbase=max(GI.MP[index,2])
  }

    x0 <- as.numeric(GI.MP[,2])
    y0 <- as.numeric(GI.MP[,3])
    z0 <- as.numeric(GI.MP[,1])
	position=order(y0,decreasing = TRUE)
    index0=GAPIT.Pruning(y0[position],numSlots=numSlots)
	index=position[index0]
	x=x0[index]
	y=y0[index]
	z=z0[index]

  pdf(paste("GAPIT.", name.of.trait,".Manhattan-Plot.Genomewise.pdf" ,sep = ""), width = 10)
  par(mar = c(5,5,5,1))

 plot(y~x,xlab=expression(Chromosome),ylab=expression(-log[10](italic(p))) ,
       cex.lab=2,col=ifelse(z%%2==0,"orangered","navy"),axes=FALSE,type = "p",pch=20,main = paste(name.of.trait,sep=" "))

 axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
  axis(2, at=1:y.lim,cex.axis=1.5,labels=1:y.lim,tick=F)
  box()
  dev.off()
} #Genomewise plot

  #print("GAPIT.Manhattan accomplished successfully!")
} #end of GAPIT.Manhattan




##############################################################################################
GAPIT.Filter.SNPs.with.Low.MAFs <- function(GI = GI, P.values = P.values, maf = maf, nobs = nobs, MAF.Filter.Rate = 0.00){

#Object: Filter out SNPs with MAFs less than a user-defined threshold. 
#Output: PWI.Filtered
#Authors: Alex Lipka, Zhiwu Zhang
# Last update: May 6, 2011 
##############################################################################################

if(is.null(GI))  PWI.Filtered=NULL
          
#Merge 
if(!is.null(GI))
{
PWI.with.MAFs <- cbind(GI,P.values, maf, nobs)
colnames(PWI.with.MAFs [c(5,6)]) <- c("MAF", "N")

#Filter out SNPs with MAF less than the threshold
PWI.Filtered <- PWI.with.MAFs[which(PWI.with.MAFs[,5] >= MAF.Filter.Rate), ]
}
#print("GAPIT.Filter.SNPs.with.Low.MAFs accomplished successfully!")

return(PWI.Filtered)
}






##############################################################################################
GAPIT.Compression.Visualization <- function(Compression = Compression, name.of.trait = name.of.trait){
#Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
#Output: Three pdfs: One of the log likelihood function, one of the genetic and error variance component,
#                    and one of the heritabilities
#Authors: Alex Lipka and Zhiwu Zhang 
# Last update: May 10, 2011 
##############################################################################################

#Graph the optimum compression 

if(length(Compression)<=6) Compression=t(as.matrix(Compression[which(Compression[,4]!="NULL"&Compression[,4]!="NaN"),]))
if(length(Compression)>6) Compression=Compression[which(Compression[,4]!="NULL"&Compression[,4]!="NaN"),]

LL=as.numeric(Compression[,4])
Compression.best=Compression[which(LL==min(LL)),] 

if(length(Compression.best)>6) Compression.best=Compression.best[1,] #Keep the first if multiple combinations
variance=as.numeric(Compression.best[5:6])
colors <- c("grey50","grey70")
labels0 <- round(variance/sum(variance) * 100, 1)
labels <- paste(labels0, "%", sep="")
LL.best0=as.numeric(Compression.best[4]  )
LL.best=floor(LL.best0*100)/100
theOptimum=paste(c(Compression.best[c(1:3)],LL.best) )

pdf(paste("GAPIT.", name.of.trait,".Optimum.pdf", sep = ""), width = 14)
par(mfrow = c(1,1), mar = c(1,1,5,5), lab = c(5,5,7))
pie(variance,  col=colors, labels=labels,angle=45)
legend(1.0, 0.5, c("Genetic varaince","Residual varaiance"), cex=1.5, 
   fill=colors)

#Display the optimum compression
text(1.5,.0, "The optimum compression", col= "red")
for(i in 1:4){
text(1.5,-.1*i, theOptimum[i], col= "red")
}
dev.off() 

#Graph compression with multiple groups
if(length(unique(Compression[,3]))>1)
{
#Create a vector of colors
color.vector.basic <- c("red","blue","black", "blueviolet","indianred","cadetblue","orange")
color.vector.addition <- setdiff(c(colors()[grep("red",colors())], colors()[grep("blue",colors())]),color.vector.basic )
color.vector.addition.mixed <- sample(color.vector.addition,max(0,((length(unique(Compression[,1])) * length(unique(Compression[,2])))-length(color.vector.basic))))  
color.vector <- c(color.vector.basic,color.vector.addition.mixed )


#Create a vector of numbers for the line dot types
line.vector <-  rep(1:(length(unique(Compression[,1])) * length(unique(Compression[,2]))))



#We want to have a total of three plots, one displaying the likelihood function, one displaying the variance components, and one displaying the
# heritability 

pdf(paste("GAPIT.", name.of.trait,".Compression.multiple.group.", ".pdf", sep = ""), width = 14)

par(mfrow = c(2,3), mar = c(5,5,1,1), lab = c(5,5,7))

# Make the likelihood function plot

k <- 1
for(i in 1:length(unique(Compression[,1]))){
  for(j in 1:length(unique(Compression[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,4])  
      plot(y~x,type="l", pch = 30, lty = line.vector[i], ylim=c(min(as.numeric(Compression[,4])),max(as.numeric(Compression[,4]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
      col = color.vector[j], xlab = "Number of Groups", ylab = "-2Log Likelihoood", )
      label = paste(c(as.character(unique(Compression[,1]))[k]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,4])  
      lines(y~x,type="l", pch = 30, lty = line.vector[i], col = color.vector[j])
      label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
   }
 }
 #Make a legend
  #legend("topright",  label, fill = color.vector) 

 

# Make the genetic variance component plots
k <- 1
for(i in 1:length(unique(Compression[,1]))){
  for(j in 1:length(unique(Compression[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,5])  
      plot(y~x,type="l", pch = 17,  lty = line.vector[i], ylim=c(min(as.numeric(Compression[,5])),max(as.numeric(Compression[,5]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
      col = color.vector[j], xlab = "Number of Groups", ylab = "Genetic Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,5])  
      lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
   }
 }
 #Make a legend
  #legend("topleft",  label, fill = color.vector) 


# Make the residual variance component plots
k <- 1
for(i in 1:length(unique(Compression[,1]))){
  for(j in 1:length(unique(Compression[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,6])  
      plot(y~x,type="l", pch = 17,  ylim=c(min(as.numeric(Compression[,6])),max(as.numeric(Compression[,6]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
      col = color.vector[j], xlab = "Number of Groups", ylab = "Residual Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,6])  
      lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
   }
 }
 #Make a legend
  #legend("topright",  label, fill = color.vector) 


#calculate total variance and h2
heritablilty.vector <- as.numeric(Compression[,5])/(as.numeric(Compression[,5]) + as.numeric(Compression[,6]))
totalVariance.vector <- as.numeric(as.numeric(Compression[,5]) + as.numeric(Compression[,6]))
Compression.h2 <- cbind(Compression, heritablilty.vector,totalVariance.vector)

# Make the total variance component plots
k <- 1
for(i in 1:length(unique(Compression.h2[,1]))){
  for(j in 1:length(unique(Compression.h2[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,8])  
      plot(y~x,type="l", pch = 17,  lty = line.vector[k], ylim=c(min(as.numeric(Compression.h2[,8])),max(as.numeric(Compression.h2[,8]))), xlim = c(min(as.numeric(Compression.h2[,3])),max(as.numeric(Compression.h2[,3]))),
      col = color.vector[1], xlab = "Number of Groups", ylab = "Total Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,8]) 
      lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
   }
 }
 #Make a legend
  #legend("topright",  label, fill = color.vector) 
  

# Make the heritability plots 
k <- 1
for(i in 1:length(unique(Compression[,1]))){
  for(j in 1:length(unique(Compression[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,7]) 

      plot(y~x,type="l", pch = 17,  lty = line.vector[k], ylim=c(min(as.numeric(Compression.h2[,7])),max(as.numeric(Compression.h2[,7]))), xlim = c(min(as.numeric(Compression.h2[,3])),max(as.numeric(Compression.h2[,3]))),
      col = color.vector[1], xlab = "Number of Groups", ylab = "Heritability", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,7])  
      lines(y~x,type="l", lty = line.vector[i], pch = 17, col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }       
   }
 }
 
 #Make a legend
  #legend("topleft",  label, fill = color.vector) 
  
legend.col= 1+floor(length(unique(Compression[,1])) * length(unique(Compression[,2]))/20)
line.style=rep(1:length(unique(Compression[,1])), each = length(unique(Compression[,2])))      
line.color=rep(1:length(unique(Compression[,2])), length(unique(Compression[,1])))



# Make labels
      plot(0~0,axes=FALSE,type="l",ylab = "",xlab = "",frame.plot=FALSE)
      legend("topleft",  label, col = color.vector[line.color], lty = line.style, ncol=legend.col,horiz=FALSE) 
   
 
dev.off()
}#end of Graph compression with multiple groups

#Graph compression with single groups
if(length(unique(Compression[,3]))==1& length(unique(Compression[,1]))*length(unique(Compression[,2]))>1)
{

#Graph the compression with only one group
pdf(paste("GAPIT.Compression.single.group.", name.of.trait, ".pdf", sep = ""), width = 14)
par(mfrow = c(2,2), mar = c(5,5,1,1), lab = c(5,5,7))

nkt=length(unique(Compression[,1]))
nca=length(unique(Compression[,2]))
kvr=rep(c(1:nkt),nca)
kvc0=rep(c(1:nca),nkt)
kvc=as.numeric(t(matrix(kvc0,nca,nkt)))
kt.name=Compression[1:nkt,1]

ca.index=((1:nca)-1)*nkt+1
ca.name=Compression[ca.index,2]

KG<- t(tapply(as.numeric(Compression[,4]), list(kvr, kvc), mean))
colnames(KG)=kt.name
barplot(as.matrix(KG),  ylab= "-2 Log Likelihood",beside=TRUE, col=rainbow(length(unique(Compression[,2]))))


KG<- t(tapply(as.numeric(Compression[,5]), list(kvr, kvc), mean))
colnames(KG)=kt.name
barplot(as.matrix(KG),  ylab= "Genetic varaince", beside=TRUE, col=rainbow(length(unique(Compression[,2]))))

KG<- t(tapply(as.numeric(Compression[,6]), list(kvr, kvc), mean))
colnames(KG)=kt.name
barplot(as.matrix(KG),  ylab= "Residual varaince", beside=TRUE, col=rainbow(length(unique(Compression[,2]))))

KG<- t(tapply(as.numeric(Compression[,5])/(as.numeric(Compression[,5])+as.numeric(Compression[,6])), list(kvr, kvc), mean))
colnames(KG)=kt.name
barplot(as.matrix(KG),  ylab= "Heritability", beside=TRUE, col=rainbow(length(unique(Compression[,2]))),ylim=c(0,1))

legend("topleft", paste(t(ca.name)), cex=0.8,bty="n", fill=rainbow(length(unique(Compression[,2]))),horiz=TRUE)
dev.off() 
} #end of Graph compression with single groups

#print("GAPIT.Compression.Visualization accomplished successfully!")

}#GAPIT.Compression.Plots ends here





##############################################################################################
GAPIT.Table <- function(final.table = final.table, name.of.trait = name.of.trait,FDR.Filter.Rate=1){
#Object: Make and export a table of summary information from GWAS
#Output: A table summarizing GWAS results
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: May 10, 2011 
##############################################################################################

#Filter SNPs by FDR
index=(final.table[,7]<=FDR.Filter.Rate)
final.table=final.table[index,]

#Export this summary table as an excel file
write.table(final.table, paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)


#print("GAPIT.Table accomplished successfully!")
  

}   #GAPIT.Table ends here



##############################################################################################
GAPIT.GS.Visualization <- function(gsBLUP = gsBLUP, BINS=BINS, name.of.trait = name.of.trait){
#Object: To build heat map to show distribution of BLUP and PEV
#Output: pdf
#Authors: Zhiwu Zhang 
# Last update: May 15, 2011 
##############################################################################################

pdf(paste("GAPIT.", name.of.trait,".GS", ".pdf", sep = ""), width = 14)
par(mfrow = c(1,1), mar = c(1,1,5,5), lab = c(5,5,7))

nBin=BINS
BLUP= gsBLUP[,5]
PEV = gsBLUP[,6]

if(BLUP[1]=="NaN"){
  warning ("It was not converged. BLUP was not created!")
}
if(BLUP[1]!="NaN")
{
  range.BLUP=max(BLUP)-min(BLUP)
  range.PEV=max(PEV)-min(PEV)
  
  interval.BLUP=range.BLUP/nBin
  interval.PEV=range.PEV/nBin
  
  
  bin.BLUP=floor(BLUP/max(BLUP)*nBin)*max(BLUP)/nBin
  bin.PEV=floor(PEV/max(PEV)*nBin)*max(PEV)/nBin
  
  distinct.BLUP=unique(bin.BLUP)
  distinct.PEV=unique(bin.PEV)
  
  Position.BLUP=match(bin.BLUP,distinct.BLUP,nomatch = 0)
  Position.PEV=match(bin.PEV,distinct.PEV,nomatch = 0)
  
  value=matrix(1,length(Position.BLUP))
  KG<- (tapply(as.numeric(value), list(Position.BLUP, Position.PEV), sum))
  
  rownames(KG)=distinct.BLUP
  colnames(KG)=distinct.PEV
  
  length(unique(bin.BLUP))
  length(unique(bin.PEV))
  
  
  nba_heatmap <- heatmap(KG, Rowv=NA, Colv=NA, 
  col = cm.colors(256), scale="column", margins=c(5,10))
  dev.off() 
}
#print("GAPIT.GS.Visualization accomplished successfully!")

}   #GAPIT.GS.Visualization ends here


      
##############################################################################################
GAPIT.kinship.loiselle <- function(snps, method="additive", use="all") {
# Object: To calculate the kinship matrix using the method of Loiselle et al. (1995)
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 


  #Number of SNP types that are 0s
  n0 <- sum(snps==0,na.rm=TRUE)
  #Number of heterozygote SNP types
  nh <- sum(snps==1,na.rm=TRUE)
  #Number of SNP types that are 1s
  n1 <- sum(snps==2,na.rm=TRUE)
  #Number of SNP types that are missing
  nNA <- sum(is.na(snps))
  

 
  #Self explanatory
  dim(snps)[1]*dim(snps)[2]
  #stopifnot(n0+nh+n1+nNA == length(snps))

    
  #Note that the two lines in if(method == "dominant") and if(method == "recessive") are found in
  #if(method == "additive").  Worry about this only if you have heterozygotes, which you do not.
  if ( method == "dominant" ) {
    flags <- matrix(as.double(colMeans(snps,na.rm=TRUE) > 1),ncol(snps),nrow(snps))
    snps[!is.na(snps) && (snps == 1)] <- flags[!is.na(snps) && (snps == 1)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(colMeans(snps,na.rm=TRUE) < 1),ncol(snps),nrow(snps))
    snps[!is.na(snps) && (snps == 1)] <- flags[!is.na(snps) && (snps == 1)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(colMeans(snps,na.rm=TRUE) > 1),ncol(snps),nrow(snps))
    dsnps[!is.na(snps) && (snps==1)] <- flags[is.na(snps) && (snps==1)]
    flags <- matrix(as.double(colMeans(snps,na.rm=TRUE) < 1),ncol(snps),nrow(snps))
    rsnps[!is.na(snps) && (snps==1)] <- flags[is.na(snps) && (snps==1)]
    snps <- cbind(dsnps,rsnps)
  }

  #mafs is a (# lines)x(# SNPs)x matrix.  The rows mafs are identical, and the ij^th element is the average
  #allele frequency for the SNP in the j^th column.
  
  #if(use == "all") imputes missing SNP type values with the expected (average) allele frequency.
  if ( use == "all" ) {
    mafs <- matrix(colMeans(snps,na.rm=TRUE),ncol(snps),nrow(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    mafs <- matrix(colMeans(snps,na.rm=TRUE),ncol(snps),nrow(snps))
    snps <- snps[colSums(is.na(snps))==0,]
  }
  mafs_comp <- 1-mafs
  snps_comp <- 1-snps
  

  n <- nrow(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1
  #Create the k term on page 1422 of Loiselle et al. (1995)

  missing <- rep(NA, dim(snps)[2])  
  for(i in 1:dim(snps)[2]) {
    missing[i] <- sum(is.na(snps[,i]))
  }
  

  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      Num_First_Term_1 <- (snps[i,]-mafs[i,])*(snps[j,]-mafs[j,])
      Num_First_Term_2 <- (snps_comp[i,]-mafs_comp[i,])*(snps_comp[j,]-mafs_comp[j,])
      First_Term <- sum(Num_First_Term_1)+sum(Num_First_Term_2)

      Num_Second_Term_1 <- mafs[,i]*(1-mafs[,i])
      Num_Second_Term_2 <- mafs_comp[,i]*(1-mafs_comp[,i])
      Num_Second_Term_Bias_Correction <- 1/((2*n)-missing - 1)
      Num_Second_Term <-  Num_Second_Term_1 + Num_Second_Term_2
      Second_Term <- sum(Num_Second_Term*Num_Second_Term_Bias_Correction)

      Third_Term <- sum(Num_Second_Term) 
      
      f <- (First_Term + Second_Term)/Third_Term
      if(f <= 0){
            K[i,j] <- 0
      }
      else{
            K[i,j] <- f
      }
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}


##############################################################################################
GAPIT.PCA <- function(X,taxa, PC.number = min(ncol(X),nrow(X))){
# Object: Conduct a principal component analysis, and output the prinicpal components into the workspace,
#         a text file of the principal components, and a pdf of the scree plot
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 

#Conduct the PCA 
PCA.X <- prcomp(X)

#Create a Scree plot 

pdf("GAPIT.Scree.Plot.pdf", width = 12, height = 12)

  screeplot(PCA.X, type="lines")

dev.off()


#Write the PCs into a text file
PCs <- cbind(taxa,as.data.frame(PCA.X$x[,1:PC.number]))

#Remove duplicate

PCs.unique <- unique(PCs[,1])
PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]


write.table(PCs, "GAPIT.Principal.Components.csv", quote = FALSE, sep = ",", row.names = TRUE,col.names = TRUE)

#Return the PCs

return(PCs)


}



      
##############################################################################################
GAPIT.kinship.loiselle <- function(snps, method="additive", use="all") {
# Object: To calculate the kinship matrix using the method of Loiselle et al. (1995)
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 
 
  #Number of SNP types that are 0s
  n0 <- sum(snps==0,na.rm=TRUE)
  #Number of heterozygote SNP types
  nh <- sum(snps==0.5,na.rm=TRUE)
  #Number of SNP types that are 1s
  n1 <- sum(snps==1,na.rm=TRUE)
  #Number of SNP types that are missing
  nNA <- sum(is.na(snps))
  

 
  #Self explanatory
  dim(snps)[1]*dim(snps)[2]
  #stopifnot(n0+nh+n1+nNA == length(snps))

    
  #Note that the two lines in if(method == "dominant") and if(method == "recessive") are found in
  #if(method == "additive").  Worry about this only if you have heterozygotes, which you do not.
  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  #mafs is a (# SNPs)x(# lines) matrix.  The columns of mafs are identical, and the ij^th element is the average
  #allele frequency for the SNP in the i^th row.
  
  #if(use == "all") imputes missing SNP type values with the expected (average) allele frequency.
  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps <- snps[rowSums(is.na(snps))==0,]
  }
  mafs_comp <- 1-mafs
  snps_comp <- 1-snps
  

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1
  #Create the k term on page 1422 of Loiselle et al. (1995)

  missing <- rep(NA, dim(snps)[1])  
  for(i in 1:dim(snps)[1]) {
    missing[i] <- sum(is.na(snps[i,]))
  }
  

  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      Num_First_Term_1 <- (snps[,i]-mafs[,i])*(snps[,j]-mafs[,j])
      Num_First_Term_2 <- (snps_comp[,i]-mafs_comp[,i])*(snps_comp[,j]-mafs_comp[,j])
      First_Term <- sum(Num_First_Term_1)+sum(Num_First_Term_2)

      Num_Second_Term_1 <- mafs[,i]*(1-mafs[,i])
      Num_Second_Term_2 <- mafs_comp[,i]*(1-mafs_comp[,i])
      Num_Second_Term_Bias_Correction <- 1/((2*n)-missing - 1)
      Num_Second_Term <-  Num_Second_Term_1 + Num_Second_Term_2
      Second_Term <- sum(Num_Second_Term*Num_Second_Term_Bias_Correction)

      Third_Term <- sum(Num_Second_Term) 
      
      f <- (First_Term + Second_Term)/Third_Term

      K[i,j] <- f
      if(K[i,j]<0) K[i,j]=0
      
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}

##############################################################################################
GAPIT.replaceNaN <- function(LL) {
#handler of grids with NaN log
#Authors: Zhiwu Zhang
# Last update: may 12, 2011 
##############################################################################################

#handler of grids with NaN log 
index=(LL=="NaN")
if(length(index)>0) theMin=min(LL[!index])
if(length(index)<1) theMin="NaN"
LL[index]=theMin
return(LL)    
}

################################################################################################################
emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
              esp=1e-10, eig.L = NULL, eig.R = NULL) {
# Authors: Hyun Min Kang
# Modified (only one line) by Zhiwu Zhang to handle non-defined LL ("NaN") by replacing it with the worst LL.
# Last update: June 8, 2011 
################################################################################################################


  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
  
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
  
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  
  #handler of grids with NaN log
  optLL=GAPIT.replaceNaN(optLL)   
  
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta

  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

##############################################################################################
GAPIT.Sampler <- function(GD=NULL,dataPath=NULL,numFiles=1,GFile=NULL,GFileExt=NULL,seed=0,ratio=1,model="Add",
                          genoFormat=NULL, GDFile=NULL, GDFileExt=NULL){
#Object: To sample genotype from multiple genotype files
#Output: genotype data sampled
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: June 11, 2011
##############################################################################################

#Take a random sample of SNPs from the first file. Output is "GD.Ssample"
set.seed(seed)
n= ncol(GD)
sample=sample(1:n,floor(n*ratio))
GD.Sample=GD[,sample]
rm(GD)
gc()
#Take a random sample of SNPs from the remaining files
if (numFiles>1)
{
  if(genoFormat=="Hapmap"){
    for (file in 2:numFiles)
    {
      G <- read.table(paste(dataPath,GFile,file, ".",GFileExt,sep=""), head = FALSE)

      hm=GAPIT.HapMap(G,model=model)
      rm(G)
      gc()
      n= ncol(hm$GD)
      set.seed(seed+file)
      sample=sample(1:n,floor(n*ratio))
      GD.Sample=cbind(GD.Sample,hm$GD[sample])
      rm(hm)
      gc()
    } # end of multiple genotype files
  } #end of numfile>1
  if(genoFormat=="EMMA"){
    for (file in 2:numFiles)
    {
      GD <- read.table(paste(dataPath,GDFile,file, ".",GDFileExt,sep=""), head = TRUE)
      GD=GD[,-1]
      n= ncol(GD)
      set.seed(seed+file)
      sample=sample(1:n,floor(n*ratio))
      GD.Sample=cbind(GD.Sample,GD[sample])
      rm(GD)
      gc()
    } # end of multiple genotype files
  } #end of numfile>1
}
return(GD.Sample)
}




