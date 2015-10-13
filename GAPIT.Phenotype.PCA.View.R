`GAPIT.Phenotype.PCA.View` <-function(PC=NULL,myY=NULL){
# Object: Analysis PCA effection for Phenotype data ,result:a pdf of the scree plot
# myG:Genotype data
# myY:Phenotype data

# Authors: You Tang
# Last update: Sep 7, 2015 
############################################################################################## 
if(is.null(PC)){stop("Validation Invalid. Please input four PC value  !")}
if(is.null(myY)){stop("Validation Invalid. Please select read valid Phenotype flies  !")}

y<-myY[!is.na(myY[,2]),c(1:2)]

traitname=colnames(y)[2]

cv1<-PC[!is.na(match(PC[,1],y[,1])),]
y1<-y[!is.na(match(y[,1],cv1[,1])),]

y2<-y1[order(y1[,1]),]
cv2<-cv1[order(cv1[,1]),]

pdf(paste("PCA_Phenotype view_",traitname,".pdf",seq=""), width =9, height = 6)
par(mar = c(5,5,5,5))
par(mfrow=c(2,2))
plot(y2[,2],cv2[,2],bg="lightgray",xlab="Phenotype",ylab="PC1",main="",cex.lab=1.4)
plot(y2[,2],cv2[,3],bg="lightgray",xlab="Phenotype",ylab="PC2",main="",cex.lab=1.4)
plot(y2[,2],cv2[,4],bg="lightgray",xlab="Phenotype",ylab="PC3",main="",cex.lab=1.4)
plot(y2[,2],cv2[,5],bg="lightgray",xlab="Phenotype",ylab="PC4",main="",cex.lab=1.4)
dev.off()


print(paste("GAPIT.Phenotype.PCA.View ", ".output pdf generate.","successfully!" ,sep = ""))

#GAPIT.Phenotype.View
}
