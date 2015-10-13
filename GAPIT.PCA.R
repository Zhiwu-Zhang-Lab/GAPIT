`GAPIT.PCA` <-
function(X,taxa, PC.number = min(ncol(X),nrow(X)),file.output=TRUE){
# Object: Conduct a principal component analysis, and output the prinicpal components into the workspace,
#         a text file of the principal components, and a pdf of the scree plot
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 

#Conduct the PCA 
print("Calling prcomp...")
PCA.X <- prcomp(X)

print("Creating PCA graphs...")
#Create a Scree plot 
if(file.output & PC.number>1) {
pdf("GAPIT.PCA.eigenValue.pdf", width = 12, height = 12)
par(mar = c(5,5,5,5))
screeplot(PCA.X, type="lines")
dev.off()

pdf("GAPIT.PCA.pdf", width = 8, height = 8)
par(mar = c(5,5,5,5))
maxPlot=min(as.numeric(PC.number[1]),3)

for(i in 1:(maxPlot-1)){
for(j in (i+1):(maxPlot)){
plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i,sep=""),ylab=paste("PC",j,sep=""),pch=19,col="red",cex.axis=1.3,cex.lab=1.4, cex.axis=1.2, lwd=2,las=1)

}
}
dev.off()

#output 3D plot
if(ncol(X)>2){

pdf("GAPIT.PCA.3D.pdf", width = 7, height = 7)
par(mar = c(5,5,5,5))
scatterplot3d(PCA.X$x[,1],PCA.X$x[,2],PCA.X$x[,3],xlab=paste("PC",1,sep=""),ylab=paste("PC",2,sep=""),zlab=paste("PC",3,sep="") ,pch=20,color="red",col.axis="blue",cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3,angle=55,scale.y=0.7)
dev.off()
}
}
print("Joining taxa...")
#Extract number of PCs needed
PCs <- cbind(taxa,as.data.frame(PCA.X$x))

#Remove duplicate (This is taken care by QC)
#PCs.unique <- unique(PCs[,1])
#PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]

eigenvalues <- PCA.X$sdev^2

print("Exporting PCs...")
#Write the PCs into a text file
if(file.output) write.table(PCs, "GAPIT.PCA.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) write.table(PCA.X$rotation[,1:PC.number], "GAPIT.PCA.loadings.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) write.table(eigenvalues, "GAPIT.PCA.eigenvalues.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Return the PCs
return(list(PCs=PCs,EV=PCA.X$sdev^2,nPCs=NULL))
}
