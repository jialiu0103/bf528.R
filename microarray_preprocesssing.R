
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")

library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")

##install and libaray packages 



normal_cels <- '/Users/liujia/Desktop/intern/datasets/NnC'
affy_norm <- ReadAffy(celfile.path=normal_cels)
View(affy_norm)
norm_arry <- rma(affy_norm)
getClass(affy_norm)
##read in array and normalize 

##Relative Log Expression, Normalized Unscaled Standard Error
##compute RLE and NUSE scores of the microarray samples
my_plm<- fitPLM(affy_norm,normalize = T, background = T)
Mbox(my_plm,main="RLE") #this is a boxplot of rle. cause i want to take a look at distribution

sta_rle <- RLE(my_plm,type="stats")
my_med <- sta_rle[1,]
#plot histgram
hist_rle <- hist(my_med)

##nuse histgram
sta_nuse <- NUSE(my_plm,type="stats")
my_med_nuse <- sta_nuse[1,]
hist_nuse <- hist(my_med_nuse)


