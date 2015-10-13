`GAPIT.Bread` <-
function(Y=NULL,CV=NULL,Z=NULL,KI=NULL,GK=NULL,GD=NULL,GM=NULL,
              method=NULL,delta=NULL,vg=NULL,ve=NULL,LD=0.01,GTindex=NULL,
              file.output=TRUE){
#Object: To calculate p-values of SNPs by using method of GLM, MLM, CMLM, FaST, SUPER and DC  
#Straitegy: NA
#Output: GWAS, GPS,REMLs,vg,ve,delta
#intput: 
#Y: phenotype with columns of taxa,Y1,Y2...
#CV: covariate variables with columns of taxa, v1,v2...
#GD: same as GK. This is the genotype to screen, the columns are taxa,SNP1,SNP2,...
#GK: Genotype data in numerical format, taxa goes to row and snp go ti columns. the first column is taxa
#GM: Genotype map with columns of snpID,chromosome and position
#method: Options are GLM, MLM, CMLM, FaST, SUPER and DC 
#Authors: Zhiwu Zhang
#Last update: November 2, 2011
##############################################################################################
#print("GAPIT.SUPER in progress...")

#Performing first screening with GLM
if(method=="GLM"){
#print("---------------screening by GLM----------------------------------")

  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=0,			
  group.to=0,
  QC=FALSE,
  GTindex=GTindex,
  file.output=file.output				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with MLM
if(method=="MLM"){
#print("---------------screening by MLM----------------------------------")

  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=nrow(Y),			
  group.to=nrow(Y),
  QC=FALSE,
  GTindex=GTindex,
  file.output=file.output				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with Compressed MLM
if(method=="CMLM"){
#print("---------------screening by CMLM----------------------------------")
  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=1,			
  group.to=nrow(Y),
  QC=FALSE,
  GTindex=GTindex,
  file.output=file.output				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with FaST-LMM
if(method=="FaST" | method=="SUPER"| method=="DC")
{
  GWAS=NULL
  GPS=NULL
  if(!is.null(vg) & !is.null(vg) & is.null(delta)) delta=ve/vg
  if(is.null(vg) & is.null(ve))
  {

    myFaSTREML=GAPIT.get.LL(pheno=matrix(Y[,-1],nrow(Y),1),geno=NULL,snp.pool=as.matrix(GK[,-1]),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])))
    
#print("Transfer data...")    
    REMLs=-2*myFaSTREML$LL  
    delta=myFaSTREML$delta
    vg=myFaSTREML$vg
    ve=myFaSTREML$ve
    #GPS=myFaSTREML$GPS
  }

mySUPERFaST=GAPIT.SUPER.FastMLM(ys=matrix(Y[,-1],nrow(Y),1),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])),snp.pool=as.matrix(GK[-1]), xs=as.matrix(GD[GTindex,-1]),vg=vg,delta=delta,LD=LD,method=method)
GWAS=cbind(GM,mySUPERFaST$ps,mySUPERFaST$stats,mySUPERFaST$dfs,mySUPERFaST$effect)
}#End of if(method=="FaST" | method=="SUPER")

#print("GAPIT.Bread succeed!")  
return (list(GWAS=GWAS, GPS=GPS,REMLs=REMLs,vg=vg,ve=ve,delta=delta))
} #end of GAPIT.Bread

