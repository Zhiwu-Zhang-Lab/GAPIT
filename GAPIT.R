`GAPIT` <-
function(Y=NULL,G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,CV.Inheritance=NULL,GP=NULL,GK=NULL,
                group.from=30 ,group.to=1000000,group.by=10,DPP=100000, 
                kinship.cluster="average", kinship.group='Mean',kinship.algorithm="VanRaden",                                                    
                bin.from=10000,bin.to=10000,bin.by=10000,inclosure.from=10,inclosure.to=10,inclosure.by=10,
                SNP.P3D=TRUE,SNP.effect="Add",SNP.impute="Middle",PCA.total=0, SNP.fraction = 1, seed = 123, BINS = 20,SNP.test=TRUE, 
                SNP.MAF=0,FDR.Rate = 1, SNP.FDR=1,SNP.permutation=FALSE,SNP.CV=NULL,SNP.robust="GLM",                             
                file.from=1, file.to=1, file.total=NULL, file.fragment = 99999,file.path=NULL, 
                file.G=NULL, file.Ext.G=NULL,file.GD=NULL, file.GM=NULL, file.Ext.GD=NULL,file.Ext.GM=NULL, 
                ngrid = 100, llim = -10, ulim = 10, esp = 1e-10,
                LD.chromosome=NULL,LD.location=NULL,LD.range=NULL,
                sangwich.top=NULL,sangwich.bottom=NULL,QC=TRUE,GTindex=NULL,LD=0.1,
                file.output=TRUE,cutOff=0.01, Model.selection = FALSE,output.numerical = FALSE,
                output.hapmap = FALSE, Create.indicator = FALSE,
				QTN=NULL, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty", Major.allele.zero = FALSE,
        method.GLM="fast.lm",method.sub="reward",method.sub.final="reward",method.bin="static",bin.size=c(1000000),bin.selection=c(10,20,50,100,200,500,1000),
        memo="",Prior=NULL,ncpus=1,maxLoop=3,threshold.output=.01,
        WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),maxOut=100,QTN.position=NULL,
        converge=1,iteration.output=FALSE,acceleration=0,iteration.method="accum",PCA.View.output=TRUE,Geno.View.output=TRUE,plot.style="Oceanic",SUPER_GD=NULL,SUPER_GS=FALSE){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Designed by Zhiwu Zhang
#Writen by Alex Lipka, Feng Tian ,You Tang and Zhiwu Zhang
#Last update: Oct 23, 2015  by Jiabo Wang add REML threshold and SUPER GK
##############################################################################################
print("--------------------- Welcome to GAPIT ----------------------------")
  
echo=TRUE
GAPIT.Version=GAPIT.0000()

Timmer=GAPIT.Timmer(Infor="GAPIT")
Memory=GAPIT.Memory(Infor="GAPIT")

#Genotype processing and calculation Kin and PC
#First call to genotype to setup genotype data

storage_PCA.total<-PCA.total
#if(PCA.total>0){
#if(PCA.total<=3){PCA.total=4}
#}

#BUS algorithm
if(kinship.algorithm=="FARM-CPU") return (GAPIT.BUS(Y=Y,GDP=GD,GM=GM,CV=CV,
  method.GLM=method.GLM,method.sub=method.sub,method.sub.final=method.sub.final,method.bin=method.bin,
  bin.size=bin.size,bin.selection=bin.selection,file.output=file.output,
  cutOff=cutOff,DPP=DPP,memo=memo,Prior=Prior,ncpus=ncpus,maxLoop=maxLoop,
  kinship.algorithm=kinship.algorithm,GP=GP,threshold.output=threshold.output,
  WS=WS,alpha=alpha,maxOut=maxOut,QTN.position=QTN.position,converge=converge,
  iteration.output=iteration.output,acceleration=acceleration,iteration.method=iteration.method))

myGenotype<-GAPIT.Genotype(G=G,GD=GD,GM=GM,KI=KI,kinship.algorithm=kinship.algorithm,PCA.total=PCA.total,SNP.fraction=SNP.fraction,SNP.test=SNP.test,
                file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G, 
                file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
                SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,
                LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,
                GP=GP,GK=GK,bin.size=NULL,inclosure.size=NULL, Timmer = Timmer,Memory=Memory,
                sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,GTindex=NULL,file.output=file.output, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero,Geno.View.output=Geno.View.output)

Timmer=myGenotype$Timmer
Memory=myGenotype$Memory

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype for all")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype for all")


KI=myGenotype$KI
PC=myGenotype$PC


genoFormat=myGenotype$genoFormat
hasGenotype=myGenotype$hasGenotype
byFile=myGenotype$byFile
fullGD=myGenotype$fullGD
GD=myGenotype$GD
GI=myGenotype$GI
GT=myGenotype$GT
G=myGenotype$G

rownames(GD)=GT
colnames(GD)=GI[,1]

if(output.numerical) write.table(GD,  "GAPIT.Genotype.Numerical.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
if(output.hapmap) write.table(myGenotype$G,  "GAPIT.Genotype.hmp.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

#In case of null Y and null GP, return genotype only  
if(is.null(Y) & is.null(GP)) return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=myGenotype$KI,PCA=myGenotype$PC,GD=data.frame(cbind(as.data.frame(GT),as.data.frame(GD))),GI=GI,G=myGenotype$G))

#In case of null Y, return genotype only          
if(is.null(Y)) return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=myGenotype$KI,PCA=myGenotype$PC,GD=data.frame(cbind(as.date.frame(GT),as.data.frame(GD))),Gi=GI,G=myGenotype$G))

rm(myGenotype)
gc()


PCA.total<-storage_PCA.total

print("--------------------Processing traits----------------------------------")
if(!is.null(Y)){
print("Phenotype provided!")
if(ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")

for (trait in 2: ncol(Y))  {
traitname=colnames(Y)[trait]

###Statistical distributions of phenotype
if(!is.null(Y) & file.output)ViewPhenotype<-GAPIT.Phenotype.View(myY=Y[,c(1,trait)],traitname=traitname,memo=memo)


###Correlation between phenotype and principal components
if(!is.null(Y)&!is.null(PC) & file.output & PCA.total>0 & PCA.View.output){

myPPV<-GAPIT.Phenotype.PCA.View(
PC=PC,
myY=Y[,c(1,trait)]
)

}

print(paste("Processing trait: ",traitname,sep=""))
if(!is.null(memo)) traitname=paste(memo,".",traitname,sep="")
gapitMain <- GAPIT.Main(Y=Y[,c(1,trait)],G=G,GD=GD,GM=GM,KI=KI,Z=Z,CV=CV,CV.Inheritance=CV.Inheritance,GP=GP,GK=GK,SNP.P3D=SNP.P3D,kinship.algorithm=kinship.algorithm,
                      bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,
				              group.from=group.from,group.to=group.to,group.by=group.by,kinship.cluster=kinship.cluster,kinship.group=kinship.group,name.of.trait=traitname,
                        file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G,file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM, 
                        SNP.MAF= SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total,GAPIT.Version=GAPIT.Version,
                        GT=GT, SNP.fraction = SNP.fraction, seed = seed, BINS = BINS,SNP.test=SNP.test,DPP=DPP, SNP.permutation=SNP.permutation,
                        LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,SNP.CV=SNP.CV,SNP.robust=SNP.robust,
                        genoFormat=genoFormat,hasGenotype=hasGenotype,byFile=byFile,fullGD=fullGD,PC=PC,GI=GI,Timmer = Timmer, Memory = Memory,
                        sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,QC=QC,GTindex=GTindex,LD=LD,file.output=file.output,cutOff=cutOff, 
                        Model.selection = Model.selection, Create.indicator = Create.indicator,
						            QTN=QTN, QTN.round=QTN.round,QTN.limit=QTN.limit, QTN.update=QTN.update, QTN.method=QTN.method, Major.allele.zero=Major.allele.zero,
                        QTN.position=QTN.position,plot.style=plot.style,SUPER_GS=SUPER_GS)  
}# end of loop on trait

if(ncol(Y>2) &file.output)
{
Timmer=gapitMain$Timmer
Memory=gapitMain$Memory

file=paste("GAPIT.", "All",".Timming.csv" ,sep = "")
write.table(Timmer, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

file=paste("GAPIT.", "All",".Memory.Stage.csv" ,sep = "")
write.table(Memory, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}

if(ncol(Y)==2) {

if (!SUPER_GS){
#Evaluate Power vs FDR and type I error
myPower=NULL
if(!is.null(gapitMain$GWAS))myPower=GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=GM,GWAS=gapitMain$GWAS)


h2= as.matrix(as.numeric(as.vector(gapitMain$Compression[,5]))/(as.numeric(as.vector(gapitMain$Compression[,5]))+as.numeric(as.vector(gapitMain$Compression[,6]))),length(gapitMain$Compression[,6]),1)
colnames(h2)=c("Heritability")
  print("GAPIT accomplished successfully for single trait. Results are saved. GWAS are returned!")
  print("It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'")

  return (list(QTN=gapitMain$QTN,GWAS=gapitMain$GWAS,h2=gapitMain$h2,Pred=gapitMain$Pred,compression=as.data.frame(cbind(gapitMain$Compression,h2)), 
  kinship.optimum=gapitMain$kinship.optimum,kinship=gapitMain$kinship,PCA=gapitMain$PC,
    FDR=myPower$FDR,Power=myPower$Power,Power.Alpha=myPower$Power.Alpha,alpha=myPower$alpha,SUPER_GD=gapitMain$SUPER_GD,P=gapitMain$P,effect.snp=gapitMain$effect.snp,effect.cv=gapitMain$effect.cv))
}else{
h2= as.matrix(as.numeric(as.vector(gapitMain$Compression[,5]))/(as.numeric(as.vector(gapitMain$Compression[,5]))+as.numeric(as.vector(gapitMain$Compression[,6]))),length(gapitMain$Compression[,6]),1)
colnames(h2)=c("Heritability")
  print("GAPIT accomplished successfully for single trait. Results are saved. GPS are returned!")
  print("It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'")

  return (list(QTN=gapitMain$QTN,GWAS=gapitMain$GWAS,h2=gapitMain$h2,Pred=gapitMain$Pred,compression=as.data.frame(cbind(gapitMain$Compression,h2)), 
  kinship.optimum=gapitMain$kinship.optimum,kinship=gapitMain$kinship,PCA=gapitMain$PC,
    SUPER_GD=gapitMain$SUPER_GD,P=gapitMain$P,effect.snp=gapitMain$effect.snp,effect.cv=gapitMain$effect.cv))

}


}else{
  print("GAPIT accomplished successfully for multiple traits. Results are saved")
  print("It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'")

  
  return (list(QTN=NULL,GWAS=NULL,h2=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=gapitMain$KI,PCA=gapitMain$PC,P=gapitMain$P,effect.snp=gapitMain$effect.snp,effect.cv=gapitMain$effect.cv))
}

}# end ofdetecting null Y
}  #end of GAPIT function
#=============================================================================================

