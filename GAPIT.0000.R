`GAPIT.0000` <-
function(){
##############################################################################################
#GAPIT: Genome Association and Prediction Integrated Tool
#Objective 1: State of art methods for high  power, accuracy and speed;
#Objective 2: User friendly by design, help documents, and web forum;
#Objective 3: Comprehensive output to interpret data and results;
#Objective 4: Informative tables and high quality figures for reports and publication;

#Methods implimented: 
# 1. GLM (Structure or Q method for GWAS, Pritchard et. al. Genetics, 2000)
# 2. MLM (Q+K, Yu et. al. Nature Genetics, 2006)
# 3. gBLUP (Marker based kinship, Zhang et. al. Journal of Animal Science, 2007)
# 4. PCA (Zhao et. al. Plos Genetics, 2007)
# 5. EMMA (Kang et. al. Genetics, 2008)
# 6. CMLM (Zhang et. al. Nature Genetics, 2010)
# 7. EMMAx (Kang et. al. Nature Genetics, 2010)
# 8. P3D (Zhang et. al. Nature Genetics, 2010)
# 9. FaST-LMM (Lippert et. al. Nature Methods, 2011)
# 10. ECMLM (Li et. al. BMC Bioogy, 2014)
# 11. SUPER (Wang et. al. PLoS One, 2014)

#Designed by Zhiwu Zhang
#Authors: Alex Lipka, Feng Tian, Qishan Wang, Xiaolei Liu, Meng Li,You Tang and Zhiwu Zhang
#Citation: Lipka AE, Tian F, Wang Q, Peiffer J, Li M, Bradbury PJ, Gore MA, Buckler ES, and Zhang Z 
#GAPIT: genome association and prediction integrated tool. Bioinformatics 2012, 28:2397-2399.
GAPIT.Version="2016.03.01 (Kinship defined by Zhiwu Zhang), http://zzlab.net/GAPIT"
return(GAPIT.Version)
}
#=============================================================================================

