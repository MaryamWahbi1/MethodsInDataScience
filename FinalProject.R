# ----------------------------------- Final Project ---------------------------------------------------


#------------------------- Start Part 1 - Initial processing of the data  -----------------------------
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


# ----------------------------------- Loading data ----------------------------------------------------
#######################################################################################################
#######################################################################################################


#setting the working directory 
setwd("~/Library/CloudStorage/OneDrive-Personal/FinalProject_207186016_318388899")

#path to the directory of the data 
input = paste(getwd(), "/", sep = "")

#loading "mat.f.coding.RData"
load(paste0(input,"mat.f.coding.RData"), verbose = "TRUE")
ls()

#loading "pheno.f.RData"
load(paste0(input,"pheno.f.RData"), verbose = "TRUE")
ls()

#loading "gene.f.RData"
load(paste0(input,"gene.f.RData"), verbose = "TRUE")
ls()


# --------------------------------- Type of organ -----------------------------------------------------
#######################################################################################################
#######################################################################################################


#the vector ts.list holds the type of organs
ts.list = as.character(unique(pheno.f$SMTSD))

#function to get raw tissue data
get.raw.tissue.edata<-function(tissue.name, mat.f.coding, pheno.f){ 
  tiss.cols.1 = which(pheno.f$SMTSD %in% tissue.name)
  mat.1 = mat.f.coding[, tiss.cols.1]
  return(mat.1)
}

#35 = index of the tissue - Liver organ
tmp.tissue = ts.list[35]
print(paste("The organ we have chosen is:" ,tmp.tissue))
print(paste0("loading ", tmp.tissue, " edata"))
reads.src1 = get.raw.tissue.edata(tmp.tissue, mat.f.coding, pheno.f)
t.reads.src = t(reads.src1)
length(t.reads.src)


# ---------------------------- Cleaning in genes data -------------------------------------------------
#######################################################################################################
#######################################################################################################


#------------- Delete genes with low values - with 80% of expression is 0.1----------------------------
#######################################################################################################
#######################################################################################################


vec.1 = apply(reads.src1 , 1, function(x) length(which( x > log(0.1+1, 2) )))
row.index = which(vec.1 > (0.8*(ncol(reads.src1 ))))
# leave just rows with expression per at least 80% of the samples
src.reads = reads.src1 [row.index, ]
length(src.reads)
deletion = length(t.reads.src) - length(src.reads)
print(paste("There was a deletion of",deletion,"genes"))
remaining = length(src.reads)
print(paste("We have", remaining, "remaining genes"))


#-------------------------- Delete genes with variance = 0 --------------------------------------------
#######################################################################################################
#######################################################################################################


var.data <- apply(src.reads, 1, var) #generate variance of each row - gene 
low.var.indxs = which(var.data == 0)
if(length(low.var.indxs) > 0)
{
  data.free = src.reads
  #now we get smaller matrix, with no genes with variance 0
  src.reads <- data.free[-low.var.indxs,] 
}
deletion_variance = remaining - length(src.reads)
#print(paste("There was a deletion of",deletion_variance,"genes"))
#print(paste("We have", length(src.reads), "remaining genes"))
#We found that the genes with low values also have low variance 


# --------------------------  Remove outliers in pheno data  ------------------------------------------
#######################################################################################################
#######################################################################################################


install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("Biobase")
BiocManager::install("GO.db")
BiocManager::install("htmlTable")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
library('WGCNA')


#------------------------------- Method 1 - removing outliers -----------------------------------------
#######################################################################################################
#######################################################################################################


remove.outliers.with.SD<-function(t.reads.src)
{
  #remove outliers
  #cluster the samples and not the genes to find outliers
  A = adjacency(t(t.reads.src), type = "distance")
  #the connectivity of each human. -1 is to remove the diagonal, the cor to itself
  k = as.numeric(apply(A,2,sum))-1
  Z.k = scale(k) #standardized k
  thresholdZ.k = -3 #standard deviation
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")#the red is the outlier
  my.outliers = which(outlierColor == "red")
  #printing the outlier samples
  my.outliers.samp = (rownames(t.reads.src))[my.outliers]
  print("outlier samples to remove")
  print(my.outliers.samp)
  #plot(my.outliers.samp)
  my.good.samples = which(outlierColor == "black")
  my.good.samples.names = (rownames(t.reads.src))[my.good.samples]
  #printing the outlier samples
  #print(my.good.samples.names)
  #this is the final mat after outliers removal
  t.reads.src = t.reads.src[my.good.samples.names, ]
  return(t.reads.src)
}
t.reads.src = remove.outliers.with.SD(t(src.reads))
tissue.edata = t(t.reads.src)
View(tissue.edata)


#------------------------------- Method 2 - removing outlier 1 -----------------------------------------
#######################################################################################################
#######################################################################################################


sampleTree = hclust(dist(t(tissue.edata)), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.7);#change this to change the size of the text
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

cut_avg <- cutree (sampleTree,k = 3)
plot(sampleTree)
rect.hclust(sampleTree , k = 3, border = 2:6)
abline(h = 91, col = 'yellow')

install.packages('dendextend', dependencies = TRUE)
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(sampleTree)
avg_col_dend <- color_branches(avg_dend_obj, h = 91)
plot(avg_col_dend)


#------------------------------- Remove outlier 2 -----------------------------------------------------
#######################################################################################################
#######################################################################################################


t.reads.src = remove.outliers.with.SD(t.reads.src)
tissue.edata = t(t.reads.src)
View(tissue.edata)
#Now we will run method 2 again

#------------------------------- Method 2 - removing outlier 2 -----------------------------------------
#######################################################################################################
#######################################################################################################


sampleTree = hclust(dist(t(tissue.edata)), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.7);#change this to change the size of the text
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

cut_avg <- cutree (sampleTree,k = 3)
plot(sampleTree)
rect.hclust(sampleTree , k = 3, border = 2:6)
abline(h = 91, col = 'yellow')

avg_dend_obj <- as.dendrogram(sampleTree)
avg_col_dend <- color_branches(avg_dend_obj, h = 91)
plot(avg_col_dend)


# --------------------------  Matrix that contains feature of samples  --------------------------------
#######################################################################################################
#######################################################################################################


#help function - return matrix that contains feature of samples
get.pheno.mat.from.samplIDs <- function(pheno.mat, samples.vec){
  indexes = which(pheno.f$SAMPID %in% samples.vec) 
  #these are the column indexes of the samples
  #check the case of more than one sample for a subject
  tissue.pheno = pheno.f[indexes, ]
  #second parameter is the order that we want
  ordering = match(tissue.pheno$SAMPID, samples.vec)
  tissue.pheno = tissue.pheno[ordering,]
  if(!identical(as.character(tissue.pheno$SAMPID), samples.vec)) {print("ERROR 2: samples
in norm.mat and pheno do not match")}
  #this is the related pheno tabel ordered accoring to the intial subject ids order, as in
  #the expression table
  return(tissue.pheno)
}
samples.ids = colnames(tissue.edata) #the ids of the subjects in this tissue
tmp.pheno = get.pheno.mat.from.samplIDs(pheno.f, samples.ids)
View(tmp.pheno)


# --------------------------  Quantile normalization   ------------------------------------------------
#######################################################################################################
#######################################################################################################


#boxplot before normalization
boxplot(tissue.edata, col = rainbow(ncol(tissue.edata)))

library(preprocessCore)

#rows are genes, columns are samples
quantile.normalize.raw.gtex <- function(edata.mat)
{
  norm_edata = normalize.quantiles(as.matrix(edata.mat))
  rownames(norm_edata) = rownames(edata.mat)
  colnames(norm_edata) = colnames(edata.mat)
  return(norm_edata)
}
tissue.edata.qn = quantile.normalize.raw.gtex(tissue.edata)
View(tissue.edata.qn) 

#boxplot after normalization
boxplot(tissue.edata.qn, col = rainbow(ncol(tissue.edata.qn)))


#--------------------------- End Part 1  - ------------------------------------------------------









#--------------------------- Start Part 2  ------------------------------------------------------


#----------------Create new data frame with features---------------------------------------------
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################


install.packages("ggplot2")
install.packages('factoextra')
library(ggplot2)
library("factoextra")
library('dplyr')



SAMPID<- tmp.pheno[2]
SMRIN<- tmp.pheno[6]
SMTSISCH<- tmp.pheno[10]
SMGEBTCH<- tmp.pheno[16]
AGE <- tmp.pheno[68]
DTHHRDY<- tmp.pheno[69]

orginal_features_df <- data.frame(SAMPID,SMRIN,SMTSISCH,SMGEBTCH,AGE,DTHHRDY)
View(orginal_features_df)


#---------------------------Cleaning data frame with features -----------------------------------------
#######################################################################################################
#######################################################################################################


features_df <- data.frame(SAMPID,SMRIN,SMTSISCH,SMGEBTCH,AGE,DTHHRDY)
features_df2  <- features_df[,-1]
rownames(features_df2) <- features_df[,1]
features_df <- features_df2
View(features_df)


#--------------------------Cleaning SMRIN column-------------------------------------------------------
#######################################################################################################
#######################################################################################################



for (x in 1:length(features_df$SMRIN)) {
  if(between(features_df$SMRIN[x],5,6.9)) features_df$SMRIN[x] = paste0("Low")
  else
    if(between(features_df$SMRIN[x],7,7.9)) features_df$SMRIN[x] = paste0("Med")
    else
      if(between(features_df$SMRIN[x],8,9.9)) features_df$SMRIN[x] = paste0("High")
      
}

features_df$SMRIN
View(features_df)


#-------------------------Cleaning SMTSISCH column-----------------------------------------------------
#######################################################################################################
#######################################################################################################


#change value of SMTSISCH column
features_df$SMTSISCH <- features_df$SMTSISCH/60
for (x in 1:length(features_df$SMTSISCH)) {
  if(between(features_df$SMTSISCH[x],0,2)) features_df$SMTSISCH[x] = paste0("0-2 Hours")
  else
    if(between(features_df$SMTSISCH[x],2,4)) features_df$SMTSISCH[x] = paste0("2-4 Hours")
    else
      if(between(features_df$SMTSISCH[x],4,8)) features_df$SMTSISCH[x] = paste0("4-8 Hours")
      else
        if(between(features_df$SMTSISCH[x],8,24)) features_df$SMTSISCH[x] = paste0("8-24 Hours")
}
View(features_df)


#------------------------------------Cleaning SMGEBTCH column------------------------------------------
#######################################################################################################
#######################################################################################################


#remove "LCSET" string in the beginning of the SMGEBTCH value
features_df$SMGEBTCH<-gsub("LCSET-","",as.character(features_df$SMGEBTCH))
features_df$SMGEBTCH <- as.numeric(features_df$SMGEBTCH)
features_df$SMGEBTCH
table(features_df$SMGEBTCH)
table_SMGEBTCH <- table(features_df$SMGEBTCH) > 3
table_SMGEBTCH
top_value_SMGEBTCH <- which(table(features_df$SMGEBTCH) > 3)
top_value_SMGEBTCH
top_value_SMGEBTCH <- as.table(top_value_SMGEBTCH)
top_value_SMGEBTCH
dimnames_SMGEBTCH <- dimnames(top_value_SMGEBTCH )
dimnames_SMGEBTCH <- unlist(dimnames_SMGEBTCH, use.names = FALSE)
dimnames_SMGEBTCH 
length(dimnames_SMGEBTCH)
top_SMGEBTCH <- rownames(features_df)[which(table(features_df$SMGEBTCH) > 3)]
top_SMGEBTCH


#------------------------------------Cleaning AGE column-----------------------------------------------
#######################################################################################################
#######################################################################################################


#replace the age value to the avg age
features_df$AGE<- substring(gsub("-","",as.character(features_df$AGE)),1,2)
features_df$AGE<- as.numeric(features_df$AGE)
features_df$AGE<- features_df$AGE+5
features_df$AGE
for (x in 1:length(features_df$AGE)) {
  if(between(features_df$AGE[x],20, 39)) features_df$AGE[x] = paste0("YOUNG Adults")
  else
    if(between(features_df$AGE[x],40, 59)) features_df$AGE[x] = paste0("Middle-aged Adults")
    else
      if(between(features_df$AGE[x],60, 69)) features_df$AGE[x] = paste0("Old Adults")
      
}

View(features_df)


#------------------------------------Cleaning DTHHRDY column-------------------------------------------
#######################################################################################################
#######################################################################################################


features_df$DTHHRDY[is.na(features_df$DTHHRDY)] <- floor(mean(features_df$DTHHRDY, na.rm = TRUE))
View(features_df)


#----------------------------------------------Genes PCA-----------------------------------------------
#######################################################################################################
#######################################################################################################



#PCA with all genes
pca <- prcomp(t(tissue.edata.qn),scale=TRUE)

#PCA plot 
plot(pca$x[,1],pca$x[,2])
pca.r <- pca$sdev^2
pca.p <- round(pca.r/sum(pca.r)*100,1)
barplot(pca.p,main="Plot",xlab="Principal Component",ylab="Percent Variation")
fviz_pca_ind(pca, geom="text")
fviz_pca_ind(pca, geom="point")


#---------------------------------------SMRIN PCA------------------------------------------------------
#######################################################################################################
#######################################################################################################


#SMRIN column factor
features_df$SMRIN<- factor(features_df$SMRIN, levels = c("Low", "Med","High"))

#without Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$SMRIN,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMRIN - Quality of the RNA") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$SMRIN,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMRIN - Quality of the RNA") +
  ggtitle("2D SMRIN PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------------SMTSISCH PCA---------------------------------------------------
#######################################################################################################
#######################################################################################################


#SMTSISCH column factor
features_df$SMTSISCH <- factor(features_df$SMTSISCH, 
                               levels = c("0-2 Hours","2-4 Hours","4-8 Hours","8-24 Hours"))

#without Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$SMTSISCH,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMTSISCH") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$SMTSISCH,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMTSISCH") +
  ggtitle("2D SMTSISCH PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------------SMGEBTCH PCA---------------------------------------------------
#######################################################################################################
#######################################################################################################



key_vector <- as.vector(unlist(dimnames(table(features_df$SMGEBTCH))))
key_vector

value_vector <- as.vector(table(features_df$SMGEBTCH))
value_vector


my_vector <- features_df$SMGEBTCH
my_vector


my_new_vector <- c()

for(x in 1:length(my_vector)){
  my_new_vector[x] = value_vector[match(my_vector[x], key_vector)]
}

my_new_vector

SMGEBTCH_occurrences <- my_new_vector



#--------------------------------------- All the samples ----------------------------------------------


#without Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=SMGEBTCH_occurrences,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Repetitions experiment") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=SMGEBTCH_occurrences,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Repetitions experiment") +
  ggtitle("2D SMTSISCH PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#--------------------------------------- top_SMGEBTCH (with > 2) ----------------------------------------------


top_SMGEBTCH_2 <- rownames(features_df)[which(SMGEBTCH_occurrences  > 2)]
top_SMGEBTCH_2


df.subset_2 <- subset(tissue.edata.qn, select = top_SMGEBTCH_2 )
View(df.subset_2)

SMGEBTCH_pca_2 <- prcomp(t(df.subset_2),scale=TRUE)

top_value_SMGEBTCH_2<- SMGEBTCH_occurrences[SMGEBTCH_occurrences > 2 ]
top_value_SMGEBTCH_2



#without Ellipses
fviz_pca_ind(SMGEBTCH_pca_2 , pointshape = 16, 
             pointsize = 2, 
             habillage=top_value_SMGEBTCH_2,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMGEBTCH") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(SMGEBTCH_pca_2, pointshape = 16, 
             pointsize = 2, 
             habillage=top_value_SMGEBTCH_2,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMGEBTCH") +
  ggtitle("2D SMTSISCH PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#--------------------------------------- top_8_SMGEBTCH (with > 3 ) ----------------------------------------------

#subset data frame with top SMGEBTCH 
df.subset <- subset(tissue.edata.qn, select = top_SMGEBTCH )
View(df.subset)

#prcomp with top_SMGEBTCH subset data frame
SMGEBTCH_pca <- prcomp(t(df.subset),scale=TRUE)

#without Ellipses
fviz_pca_ind(SMGEBTCH_pca, pointshape = 16, 
             pointsize = 2, 
             habillage=dimnames_SMGEBTCH,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "SMGEBTCH") +
  ggtitle("2D SMGEBTCH PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------------AGE PCA--------------------------------------------------------
#######################################################################################################
#######################################################################################################

#AGE column factor
features_df$AGE <- factor(features_df$AGE, levels = c("YOUNG Adults", "Middle-aged Adults","Old Adults"))

#without Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$AGE,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "AGE") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$AGE,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "AGE") +
  ggtitle("2D AGE PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------------DTHHRDY PCA----------------------------------------------------
#######################################################################################################
#######################################################################################################


#without Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$DTHHRDY,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "DTHHRDY") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#add Ellipses
fviz_pca_ind(pca, pointshape = 16, 
             pointsize = 2, 
             habillage=features_df$DTHHRDY,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "DTHHRDY") +
  ggtitle("2D DTHHRDY PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------Relation between 5 features -------------------------------
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

#------------------------------Create new feature_relation data frame---------------------------
################################################################################################
################################################################################################


features_relation_df <- orginal_features_df
View(features_relation_df)


#--------------------------------------Cleaning data frame with features -----------------------
################################################################################################
################################################################################################


features_df2  <- features_relation_df[,-1]
rownames(features_df2) <- features_relation_df[,1]
features_relation_df <- features_df2
View(features_relation_df)


#------------------------------------Cleaning SMGEBTCH column-----------------------------------
################################################################################################
################################################################################################

#remove "LCSET" string in the beginning of the SMGEBTCH value
features_relation_df$SMGEBTCH <- SMGEBTCH_occurrences
View(features_relation_df)


#------------------------------------Cleaning AGE column----------------------------------------
################################################################################################
################################################################################################


#replace the age value to the avg age
features_relation_df$AGE<- substring(gsub("-","",as.character(features_relation_df$AGE)),1,2)
features_relation_df$AGE<- as.numeric(features_relation_df$AGE)
features_relation_df$AGE<- features_relation_df$AGE+5
features_relation_df$AGE
View(features_relation_df)


#------------------------------------Cleaning DTHHRDY column------------------------------------
################################################################################################
################################################################################################


features_relation_df$DTHHRDY[is.na(features_relation_df$DTHHRDY)] <- floor(mean(features_relation_df$DTHHRDY, na.rm = TRUE))
View(features_relation_df)


#------------------------------Correlation Matrix between features -----------------------------
################################################################################################
################################################################################################

#---------------------------------------Method 1------------------------------------------------

mydata.cor = cor(features_relation_df)
summary(features_relation_df)
View(features_relation_df)
mydata.cor

#---------------------------------------Method 2------------------------------------------------


library(corrplot)
corrplot(cor(features_relation_df),method="circle")
corrplot(cor(features_relation_df),method="number")


#---------------------------------------Method 3------------------------------------------------


install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot(cor(features_relation_df))


#-------------------------------Relation between SMRIN & all features---------------------------
################################################################################################
################################################################################################



#---------------------------------Relation between SMTSISCH & SMRIN-----------------------------


library(ggpubr)

lm.model<-lm(SMTSISCH~SMRIN,data=features_relation_df)
summary(lm.model)
predict (lm.model, data.frame (SMRIN=c(100)),interval="confidence")
predict (lm.model, data.frame (SMRIN=c(100)),interval="prediction")

par(mfrow=c(2,2))
plot(lm.model)
par(mfrow=c(1,1))

smrin_vs_SMTSISCH <- ggplot(features_relation_df, aes(x=SMTSISCH, y=SMRIN))+geom_point(color="red")
smrin_vs_SMTSISCH
smrin_vs_SMTSISCH <-  smrin_vs_SMTSISCH +  geom_smooth(se = FALSE, method = lm,formula = y ~ x)
smrin_vs_SMTSISCH
smrin_vs_SMTSISCH <- smrin_vs_SMTSISCH +stat_regline_equation(label.x = 3, label.y = 7)
smrin_vs_SMTSISCH
smrin_vs_SMTSISCH <-  smrin_vs_SMTSISCH+ ggtitle("SMTSISCH vs SMRIN") +xlab("SMTSISCH") + ylab("SMRIN")
smrin_vs_SMTSISCH
smrin_vs_SMTSISCH + theme(
  plot.title = element_text(color="red", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#---------------------------------relation between SMGEBTCH & SMRIN-----------------------------


table_SMGEBTCH_2 <- table(features_relation_df$SMGEBTCH)
table_SMGEBTCH_2

my_data <- as.data.frame(table_SMGEBTCH_2)                 
my_data 

Freq_SMGEBTCH<- ggplot(my_data,aes(x = Var1,y = Freq)) + geom_bar(stat="identity", fill="#FF9999", colour="black")
Freq_SMGEBTCH <- Freq_SMGEBTCH+ theme(text = element_text(size = 8))
Freq_SMGEBTCH <-  Freq_SMGEBTCH + ggtitle("Frequency of SMGEBTCH") +xlab("SMGEBTCH") + ylab("Frequency")
Freq_SMGEBTCH
Freq_SMGEBTCH + theme(
  plot.title = element_text(color="#FF9999", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


smrin_vs_smgebtch <- ggplot(features_relation_df, aes(factor(SMGEBTCH), SMRIN)) + geom_boxplot(fill="green",alpha=0.2)
smrin_vs_smgebtch <- smrin_vs_smgebtch + ggtitle("SMGEBTCH vs SMRIN") +xlab("SMGEBTCH") + ylab("SMRIN")
smrin_vs_smgebtch + theme(
  plot.title = element_text(color="green", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)

#---------------------------------relation between AGE & SMRIN----------------------------------


smrin_vs_age <- ggplot(features_relation_df, aes(factor(AGE), SMRIN)) + geom_boxplot(color="red",fill="blue",alpha=0.2)
smrin_vs_age <- smrin_vs_age + ggtitle("AGE vs SMRIN") +xlab("AGE") + ylab("SMRIN")
smrin_vs_age + theme(
  plot.title = element_text(color="red", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#---------------------------------relation between DTHHRDY & SMRIN------------------------------


smrin_vs_DTHHRDY<- ggplot(features_relation_df, aes(factor(DTHHRDY), SMRIN)) + geom_boxplot(fill="violet")
smrin_vs_DTHHRDY <- smrin_vs_DTHHRDY + ggtitle("DTHHRDY vs SMRIN") +xlab("DTHHRDY") + ylab("SMRIN")
smrin_vs_DTHHRDY + theme(
  plot.title = element_text(color="violet", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#-------------------------------Relation between SMTSISCH & all features------------------------
################################################################################################
################################################################################################


#---------------------------------relation between SMTSISCH & SMGEBTCH -------------------------



plot <- ggplot(features_relation_df, aes(factor(SMGEBTCH), SMTSISCH)) + geom_boxplot(fill="violet")
plot<- plot + ggtitle("SMGEBTCH vs SMTSISCH") +xlab("SMGEBTCH") + ylab("SMTSISCH")
plot + theme(
  plot.title = element_text(color="violet", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#---------------------------------relation between SMTSISCH & AGE ------------------------------


SMTSISCH_vs_age <- ggplot(features_relation_df, aes(factor(AGE), SMTSISCH)) + geom_boxplot(fill="purple")
SMTSISCH_vs_age <- SMTSISCH_vs_age + ggtitle("AGE vs SMTSISCH") +xlab("AGE") + ylab("SMTSISCH")
SMTSISCH_vs_age + theme(
  plot.title = element_text(color="purple", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#---------------------------------relation between SMTSISCH & DTHHRDY --------------------------


SMTSISCH_vs_DTHHRDY<- ggplot(features_relation_df, aes(factor(DTHHRDY),SMTSISCH)) + geom_boxplot(fill="#00BF7D")
SMTSISCH_vs_DTHHRDY <- SMTSISCH_vs_DTHHRDY + ggtitle("DTHHRDY vs SMTSISCH") +xlab("DTHHRDY") + ylab("SMTSISCH")
SMTSISCH_vs_DTHHRDY
SMTSISCH_vs_DTHHRDY + theme(
  plot.title = element_text(color="#00BF7D", size=14, face="bold.italic",hjust = 0.5),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold")
)


#-------------------------------Relation between SMGEBTCH & all features------------------------
################################################################################################
################################################################################################


#---------------------------------relation between SMGEBTCH & AGE ------------------------------


t <- table(features_relation_df$AGE ,features_relation_df$SMGEBTCH)

barplot(t,legend=F,beside=F,main='SMGEBTCH vs AGE',col=c(2:6))
legend("topright",
       legend = rownames(t),
       pch = 15,
       col = 2:6)


#---------------------------------relation between SMGEBTCH & DTHHRDY --------------------------

t <- table(features_relation_df$DTHHRDY ,features_relation_df$SMGEBTCH)
t

#Plot the Data
g <- ggplot(as.data.frame(t), aes(Var1, Var2)) + geom_point(aes(size = Freq), colour = "green") + theme_bw() + xlab("") + ylab("")
g <- g + ggtitle("DTHHRDY vs SMGEBTCH") +xlab("DTHHRDY") + ylab("SMGEBTCH")
g + scale_size_continuous(range=c(10,30)) + geom_text(aes(label = Freq))


#-------------------------------Relation between AGE & all features-----------------------------
################################################################################################
################################################################################################


#-------------------------------Relation between AGE & DTHHRDY ---------------------------------


t <- table(features_relation_df$AGE ,features_relation_df$DTHHRDY)

barplot(t,legend=F,beside=T,main='DTHHRDY & AGE',col=c(2:6))
legend("topright",
       legend = rownames(t),
       pch = 15,
       col = 2:6)

#------------------------------ Relation between genes_avg & all features-----------------------
################################################################################################
################################################################################################


#------------------------------ Adding avg genes in data frame-------------------------------------
################################################################################################
################################################################################################


genes_avg_df <- features_relation_df
View(genes_avg_df)

genes_avg_vector = c()
for (i in 1:117){
  genes_avg_vector[i] <-  mean(tissue.edata[,i])
}
genes_avg_vector
genes_avg_df['genes_avg'] <-genes_avg_vector
View(genes_avg_df)


#----------------------------- Correlation between genes avg and features-----------------------
################################################################################################
################################################################################################


library(corrplot)
corrplot(cor(genes_avg_df),method="circle")
corrplot(cor(genes_avg_df),method="number")


#------------------------------- SMTSISCH - Machine learning ----------------------------------
################################################################################################
################################################################################################


#---------------------------------- Random Forest - SMTSISCH ---------------------------------------------
################################################################################################
################################################################################################


install.packages("e1071")
library(e1071)


tissue.edata.cluster <- as.data.frame(t(tissue.edata.qn))
tissue.edata.cluster$cluster_by_feature <- unlist(as.factor(features_df$SMTSISCH))
View(tissue.edata.cluster)




install.packages('cowplot')
library(cowplot)
install.packages('randomForest')
library(randomForest)


tissue.edata.cluster$cluster_by_feature <- as.factor(tissue.edata.cluster$cluster_by_feature) 

model <- randomForest(cluster_by_feature ~ ., data= tissue.edata.cluster, proximity=TRUE)
model 

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=4),
  Type=rep(c("0-2 Hours","2-4 Hours","4-8 Hours","8-24 Hours"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"0-2 Hours"], 
          model$err.rate[,"2-4 Hours"], 
          model$err.rate[,"4-8 Hours"],
          model$err.rate[,"8-24 Hours"]))
ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))


model <- randomForest(cluster_by_feature~ ., data=tissue.edata.cluster, ntree=1000, proximity=TRUE)
model
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=4),
  Type=rep(c("0-2 Hours","2-4 Hours","4-8 Hours","8-24 Hours"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"0-2 Hours"], 
          model$err.rate[,"2-4 Hours"], 
          model$err.rate[,"4-8 Hours"],
          model$err.rate[,"8-24 Hours"]))
ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))


oob.values <- vector(length=10)
for(i in 1:10) {
  temp.model <- randomForest(cluster_by_feature ~ ., data=tissue.edata.cluster, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values


distance.matrix <- as.dist(1-model$proximity)
genes.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
genes.var.per <- round(genes.stuff$eig/sum(genes.stuff$eig)*100, 1)
genes.values <- genes.stuff$points
genes.data <- data.frame(Sample=rownames(genes.values),
                       X=genes.values[,1],
                       Y=genes.values[,2],
                       Status=tissue.edata.cluster$cluster_by_feature)


ggplot(data=genes.data, aes(x=X, y=Y, label=Sample)) + 
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste("Gene expression - ", genes.var.per[1], "%", sep="")) +
  ylab(paste("Gene expression  - ", genes.var.per[2], "%", sep="")) +
  ggtitle("Gene expression SMTSISCH plot using (1 - Random Forest Proximities)")


ggplot(data=genes.data, aes(x=X, y=Y, label=Sample)) + 
  geom_point(aes(color=Status),size = 4) +
  theme_bw() +
  xlab(paste("Gene expression - ", genes.var.per[1], "%", sep="")) +
  ylab(paste("Gene expression - ", genes.var.per[2], "%", sep="")) +
  ggtitle("Gene expression SMTSISCH plot using (1 - Random Forest Proximities)")

#------------------------- End Part 2 -----------------------------------------------------------------

