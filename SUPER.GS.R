`GAPIT.SUPER.GS`<-
function(Y=Y[,c(1,trait)],G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,GK=GK,kinship.algorithm=kinship.algorithm,
                      bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,
				        group.from=group.from,group.to=group.to,group.by=group.by,kinship.cluster=kinship.cluster,kinship.group=kinship.group,name.of.trait=traitname,
                        file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G,file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM, 
                        SNP.MAF= SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total,GAPIT.Version=GAPIT.Version,
                        GT=GT, SNP.fraction = SNP.fraction, seed = seed, BINS = BINS,SNP.test=SNP.test,DPP=DPP, SNP.permutation=SNP.permutation,
                        LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,SNP.CV=SNP.CV,SNP.robust=SNP.robust,
                        genoFormat=genoFormat,hasGenotype=hasGenotype,byFile=byFile,fullGD=fullGD,PC=PC,GI=GI,Timmer = Timmer, Memory = Memory,
                        sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,QC=QC,GTindex=GTindex,LD=LD,file.output=file.output,cutOff=cutOff
                        ){
#Object: To perform GPS with SUPER and Compress method
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novber 6, 2015 		
######################################################
print("--------------------- Welcome to GAPIT SUPER GS----------------------------")
Timmer=GAPIT.Timmer(Infor="GAPIT.SUPER.GS")
Memory=GAPIT.Memory(Infor="GAPIT.SUPER.GS")
my_allCV=CV
my_allKI=KI
shortcut=FALSE
LL.save=1e10
#In case of null Y and null GP, return genotype only  
thisY=Y[,2]
thisY=thisY[!is.na(thisY)]
if(length(thisY) <3){
 shortcut=TRUE
 }else{
  if(var(thisY) ==0) shortcut=TRUE
}
if(shortcut){
print(paste("Y is empty. No GWAS/GS performed for ",name.of.trait,sep=""))
return (list(compression=NULL,kinship.optimum=NULL, kinship=KI,PC=PC,GWAS=NULL, GPS=NULL,Pred=NULL, REMLs=NULL,Timmer=Timmer,Memory=Memory))
}
print("------------Examining data (QC)------------------------------------------")
if(is.null(Y)) stop ("GAPIT says: Phenotypes must exist.")
if(is.null(KI)&missing(GD) & kinship.algorithm!="SUPER") stop ("GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created.")
if(is.null(GD) & is.null(GT)) {
	GT=as.matrix(Y[,1])
	GD=matrix(1,nrow(Y),1)	
  GI=as.data.frame(matrix(0,1,3) )
  colnames(GI)=c("SNP","Chromosome","Position")
}
#merge CV with PC
if(PCA.total>0&!is.null(CV))CV=GAPIT.CVMergePC(CV,PC)
if(PCA.total>0&is.null(CV))CV=PC
if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & is.null(Z)){
taxa=as.character(Y[,1])
Z=as.data.frame(diag(1,nrow(Y)))
Z=rbind(taxa,Z)
taxa=c('Taxa',as.character(taxa))
Z=cbind(taxa,Z)
}
if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & !is.null(Z))
{
  if(nrow(Z)-1<nrow(Y)) Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)
}
noCV=FALSE
if(is.null(CV)){
noCV=TRUE
CV=Y[,1:2]
CV[,2]=1
colnames(CV)=c("taxa","overall")
}

#Remove duplicat and integragation of data
print("QC is in process...")

CVI <- CV
if(QC)
{
  qc <- GAPIT.QC(Y=Y,KI=KI, GT=GT,CV=CV,Z=Z,GK=GK)
  GTindex=qc$GTindex
  Y=qc$Y
  KI=qc$KI
  CV=qc$CV
  Z=qc$Z
  GK=qc$GK

  if(noCV)CVI=qc$CV
}
print("The value of QC is")
print(QC)
rm(qc)
gc()
#if(!is.null(sangwich.bottom)) byPass=((sangwich.bottom=="FaST" | sangwich.bottom=="SUPER" | sangwich.bottom=="DC" )& is.null(GP)   )
#if(!is.null(sangwich.top)) byPass.top=((sangwich.top=="FaST" | sangwich.top=="SUPER" | sangwich.top=="DC" )                 )

print("------------Examining data (QC) done-------------------------------------")
super_pass=FALSE
SUPER_myKI=NULL
SUPER_optimum_GD=NULL
if (!is.null(sangwich.top)) super_pass=TRUE
if(super_pass)
{
print("-------------------start SUPER BREAD-----------------------------------")
#Create GK if not provided
  if(is.null(GK)){
    nY=floor(nrow(Y)*.9)
    nG=ncol(GD)
    if(nG>nY){snpsam=sample(1:nG,nY)}else{snpsam=1:nG}
    GK=GD[GTindex,snpsam]
    SNPVar=apply(as.matrix(GK),2,var)
	#print(snpsam)
if (snpsam==1)stop ("GAPIT says: SUPER_GS must putin GD and GM.")
    GK=GK[,SNPVar>0]
    GK=cbind(as.data.frame(GT[GTindex]),as.data.frame(GK)) #add taxa
    
  }
  
  #myGD=cbind(as.data.frame(GT),as.data.frame(GD)) 
  GP=GAPIT.Bread(Y=Y,CV=CV,Z=Z,KI=KI,GK=GK,GD=cbind(as.data.frame(GT),as.data.frame(GD)),GM=GI,method=sangwich.top,GTindex=GTindex,LD=LD,file.output=file.output)$GWAS
  GK=NULL
#if(inclosure.to>nrow(Y))   ##########removed by Jiabo Wang ,unlimited number of inclosures
#{
#inclosure.to=nrow(Y)-1
#}
bin.level=seq(bin.from,bin.to,by=bin.by)
inclosure=seq(inclosure.from,inclosure.to,by=inclosure.by)
e=1 #################################number of bins and inclosure
count=0
#for (bin in bin.level){bin=bin.level[e]}
#for (inc in inclosure){inc=inclosure[e]}
for (bin in bin.level){
for (inc in inclosure){
count=count+1

  myGenotype<-GAPIT.Genotype(G=NULL,GD=NULL,GM=GI,KI=NULL,kinship.algorithm="SUPER",PCA.total=0,SNP.fraction=SNP.fraction,SNP.test=TRUE,
                    file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G, 
                    file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
                    SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,
                    LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,
                    GP=GP,GK=NULL,bin.size=bin,inclosure.size=inc,SNP.CV=SNP.CV,GTindex=GTindex,sangwich.top=NULL,sangwich.bottom=sangwich.bottom,
                    Timmer = Timmer, Memory = Memory,file.output=file.output, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)

print(paste("bin---",bin,"---inc---",inc,sep=""))
  GK=GD[GTindex,myGenotype$SNP.QTN]
  SUPER_GD=GD[,myGenotype$SNP.QTN]
  SNPVar=apply(as.matrix(GK),2,var)
  GK=GK[,SNPVar>0]
  SUPER_GD=SUPER_GD[,SNPVar>0]
  GK=cbind(as.data.frame(GT[GTindex]),as.data.frame(GK)) #add taxa
  SUPER_GD=cbind(as.data.frame(GT),as.data.frame(SUPER_GD)) #add taxa
  
  
  myBurger=GAPIT.Burger(Y=Y,CV=CV,GK=GK)  #modifed by Jiabo Wang
  myREML=myBurger$REMLs
  myVG=myBurger$vg
  myVE=myBurger$ve
print(myREML)
  
  if(count==1){
  GK.save=GK
  LL.save=myREML
  	SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD

}else{
  if(myREML<LL.save){
    GK.save=GK
    LL.save=myREML
	SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD
  }
}

  }# bin end
  }# inc end
  
  ########################BUILD SUPER KINSHIP
  ##########################################################
colnames(SUPER_optimum_GD)=c("taxa",colnames(SUPER_optimum_GD)[-1])
SUPER_taxa=as.character(SUPER_optimum_GD[,1])
SUPER_X=SUPER_optimum_GD[,-1]
SUPER_myKI_test=GAPIT.kinship.VanRaden(snps=as.matrix(SUPER_optimum_GD[,-1]))     #  build kinship
colnames(SUPER_myKI_test)=SUPER_taxa
SUPER_myKI=cbind(SUPER_taxa,as.data.frame(SUPER_myKI_test))
print("select optimum number of marker effect in GD")
print(dim(SUPER_optimum_GD))
######################################GOIN TO NEW CBLUP
Z=NULL
if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & is.null(Z)){
taxa=as.character(Y[,1])

Z=as.data.frame(diag(1,nrow(Y)))
Z=rbind(taxa,Z)
taxa=c('Taxa',as.character(taxa))
Z=cbind(taxa,Z)
}

if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & !is.null(Z))
{
  if(nrow(Z)-1<nrow(Y)) Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)

}
noCV=FALSE
if(is.null(CV)){
noCV=TRUE
CV=Y[,1:2]
CV[,2]=1
colnames(CV)=c("taxa","overall")
}
print("QC is in process...")
GK=NULL
CVI <- CV
if(QC)
{
  qc <- GAPIT.QC(Y=Y,KI=SUPER_myKI, GT=GT,CV=CV,Z=Z,GK=GK)
  GTindex=qc$GTindex
  Y=qc$Y
  KI=qc$KI
  CV=qc$CV
  Z=qc$Z
  GK=qc$GK

  if(noCV)CVI=qc$CV
}
rm(qc)
gc()
}# super_pass end

nk=1000000000
if(!is.null(KI)) nk=min(nk,nrow(KI))
if(!is.null(GK)) nk=min(nk,nrow(GK))
if(!is.null(KI))
{
  if(group.to>nk) {
    #group.to=min(nrow(KI),length(GTindex)) #maximum of group is number of rows in KI
    group.to=nk #maximum of group is number of rows in KI
    warning("The upper bound of groups is too high. It was set to the size of kinship!") 
  }
	if(group.from>nk){ 
    group.from=nk
    warning("The lower bound of groups is too high. It was set to the size of kinship!") 
  } 
}

if(!is.null(CV)){
 	if(group.to<=ncol(CV)+1) {
	#The minimum of group is number of columns in CV
	  group.from=ncol(CV)+2
	  group.to=ncol(CV)+2
	  warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!")
	}
}

  GROUP=seq(group.to,group.from,by=-group.by)#The reverse order is to make sure to include full model
if(missing("kinship.cluster")) kinship.cluster=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
if(missing("kinship.group")) kinship.group=c("Mean", "Max", "Min", "Median")
numSetting=length(GROUP)*length(kinship.cluster)*length(kinship.group)
ys=as.matrix(Y[2])
X0=as.matrix(CV[,-1])
CV.taxa=CVI[,1]
if(min(X0[,1])!=max(X0[,1])) X0 <- cbind(1, X0) #do not add overall mean if X0 has it already at first column
hold_Z=Z

 # library("EMMREML")
order_count=0
storage_reml=NULL
Compression=matrix(,numSetting,6)
colnames(Compression)=c("Type","Cluster","Group","REML","VA","VE")

for (group in GROUP)
{
  for (ca in kinship.cluster)
  {
  for (kt in kinship.group)
  {
#if(!optOnly) {print("Compressing and Genome screening..." )}
order_count=order_count+1
   #kinship.cluster="average"
   #ca=kinship.cluster[1]# need replication
   #NGROUP=2
   #group=GROUP[NGROUP] 
   #kt=kinship.group[1]

if(order_count==1)print("-------Mixed model with Kinship-----------------------------")
if(group<ncol(X0)+1) group=1 # the emma function (emma.delta.REML.dLL.w.Z) does not allow K has dim less then CV. turn to GLM (group=1)
cp <- GAPIT.Compress(KI=KI,kinship.cluster=ca,kinship.group=kt,GN=group,Timmer=Timmer,Memory=Memory)
bk <- GAPIT.Block(Z=hold_Z,GA=cp$GA,KG=cp$KG)
zc <- GAPIT.ZmatrixCompress(Z=hold_Z,GAU =bk$GA)
zrow=nrow(zc$Z)
zcol=ncol(zc$Z)-1
K = as.matrix(bk$KW)
#if (nrow(as.matrix(bk$KW))==1)
Z=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)
if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if(is.null(X0)) X0 <- matrix(1, ncol(ys), 1)

#handler of special Z and K
if(!is.null(Z)){ if(ncol(Z) == nrow(Z)) Z = NULL }
#if(!is.null(K)) {if(length(K)<2) K = NULL}
if(!is.null(K)) {if(length(K)<= 1) K = NULL}

X <-  X0 #covariate variables such as population structure
j=1
  if (is.null(Z)) Z=diag(x=1,nrow(K),ncol(K))
  if (group==1)   K=1
   emma_test <- emmreml(ys[j,], X=as.matrix(X), K=as.matrix(K), Z=Z,varbetahat=FALSE,varuhat=FALSE, PEVuhat=FALSE, test=FALSE)  
   
   print(paste(order_count, "of",numSetting,"--","Vg=",round(emma_test$Vu,4), "VE=",round(emma_test$Ve,4),"-2LL=",round(-2*emma_test$loglik,2), "  Clustering=",ca,"  Group number=", group ,"  Group kinship=",kt,sep = " "))
  emma_test_reml=-2*emma_test$loglik
  storage_reml=append(storage_reml,-2*emma_test$loglik)
Compression[order_count,1]=kt
Compression[order_count,2]=ca
Compression[order_count,3]=group
Compression[order_count,4]=-2*emma_test$loglik
Compression[order_count,5]=emma_test$Vu
Compression[order_count,6]=emma_test$Ve
  if(order_count==1){
   save_remle=emma_test_reml
   optimum_group=group
   optimum_Clustering=ca
   optimum_groupK=kt
}else{
  if(emma_test_reml<save_remle){
   save_remle=emma_test_reml
   optimum_group=group
   optimum_Clustering=ca
   optimum_groupK=kt
  }
}
}   # kt end
  } # ka end
  } # group end
 
  
cp <- GAPIT.Compress(KI=KI,kinship.cluster=optimum_Clustering,kinship.group=optimum_groupK,GN=optimum_group,Timmer=Timmer,Memory=Memory)
bk <- GAPIT.Block(Z=hold_Z,GA=cp$GA,KG=cp$KG)
zc <- GAPIT.ZmatrixCompress(Z=hold_Z,GAU =bk$GA)
zrow=nrow(zc$Z)
zcol=ncol(zc$Z)-1
K = as.matrix(bk$KW)
Z=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)
if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if(is.null(X0)) X0 <- matrix(1, ncol(ys), 1)
X <-  X0 #covariate variables such as population structure
j=1
  if (is.null(Z)) Z=diag(x=1,nrow(K),ncol(K))
  if (is.null(my_allCV)){my_allX=matrix(1,nrow(Z),1)
  }else{
    my_allX=cbind(1,as.matrix(my_allCV[,-1]))
	}
   emma_REMLE <- emmreml(ys[j,], X=as.matrix(X), K=as.matrix(K), Z=Z,varbetahat=TRUE,varuhat=TRUE, PEVuhat=TRUE, test=TRUE)  
   emma_BLUE=as.matrix(my_allX)%*%as.matrix(emma_REMLE$betahat)
   emma_BLUE=as.data.frame(cbind(as.character(my_allKI[,1]),emma_BLUE))
   colnames(emma_BLUE)=c("Taxa","emma_BLUE")
gs <- GAPIT.GS(KW=bk$KW,KO=bk$KO,KWO=bk$KWO,GAU=bk$GAU,UW=cbind(emma_REMLE$uhat,emma_REMLE$PEVuhat))
 BB= merge(gs$BLUP, emma_BLUE, by.x = "Taxa", by.y = "Taxa")
prediction=as.matrix(BB$BLUP)+as.numeric(as.vector(BB$emma_BLUE))
all_gs=cbind(BB,prediction)
  print("GAPIT SUPER GS completed successfully for multiple traits. Results are saved")
  return (list(GPS=BB,Pred=all_gs,Compression=Compression,kinship=my_allKI,SUPER_kinship=SUPER_myKI,SUPER_GD=SUPER_optimum_GD ,PC=my_allCV,Timmer=Timmer,Memory=Memory,GWAS=NULL ))

}
