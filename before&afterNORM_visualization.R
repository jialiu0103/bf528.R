library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")
library(ggplot2)
library(grid)
library(gridExtra)
##libaray packages 


##readin dataset
normal_cels='/projectnb/bf528/users/group4/others/stemaway_17dataset'
affy_norm = ReadAffy(celfile.path=normal_cels)


##namelist
namelist=colnames(affy_matrix)
namelist=''
for (i in coln){namelist=append(namelist,substring(i,22,27))}
namelist=namelist[2:35]
namelist

#visualization before normalization
affy_matrix=exprs(affy_norm)
#colnames(affy_matrix) <- factor(c(rep("control",17), rep("cancer",17)))
ta=t(affy_matrix)
ta=log(ta)
df_pca <- prcomp(ta,center = F, scale. = F)
#plot(df_pca$x[,1], df_pca$x[,2])
df_out <- as.data.frame(df_pca$x)
group <- factor(c(rep("control",17), rep("cancer",17)))
#head(df_out)

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,label = namelist))
p<-p+geom_point()+ geom_text_repel()+ ggtitle("PCA before normalization")

p
# group_list <- factor(c(rep("cancer",17), rep("control",18)))
# design <- model.matrix(~group_list)
# colnames(affy_matrix) <- levels(group_list)
# rownames(design) <- colnames(affy_matrix)

#View(affy_norm)
norm_arry <- rma(affy_norm)

#visulize after normalized
affy_matrix=exprs(norm_arry)
ta=t(affy_matrix)
df_pca <- prcomp(ta)
#plot(df_pca$x[,1], df_pca$x[,2])
df_out <- as.data.frame(df_pca$x)
group <- factor(c(rep("control",17), rep("cancer",17)))
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,label = namelist))
p<-p+geom_point()+ geom_text_repel()+ ggtitle("PCA after normalized by rma")
p

##use limma to analyse gene expression differences
library(limma)
#library(edgeR)
affy_matrix=exprs(norm_arry)
#group <- factor(c(rep("cancer",17), rep("control",18)))
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design) <- colnames(affy_matrix)

fit <- lmFit(affy_matrix, design)
fit <- eBayes(fit)
output <- topTable(fit, coef=2,n=Inf)
sum(output$adj.P.Val<0.05)
sigout=output[output$adj.P.Val<0.05,]


# I specify coef=2 because we are interested in the difference between groups, not the intercept.


## heatmap
library(pheatmap)
#affy_matrix=log(affy_matrix)
coln=colnames(affy_matrix)
dismat <- 1-cor(affy_matrix)
nm=factor(c(rep("control",17), rep("cancer",17)))
annotation_col = data.frame(SampleType = nm)
rownames(annotation_col) = namelist
colnames(dismat) = namelist
pheatmap(dismat, annotation_col = annotation_col,show_rownames=F,main = 'hierarchical clustering heatmap before normalization')



##violinplot
library(ggplot2)
library(reshape2)
mydata<-melt(affy_matrix)
p <- ggplot(mydata, aes(x=Var2, y=value))+
  ylab(expression("value, Log"[10]*"")) + 
  geom_violin()
p

