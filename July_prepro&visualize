#read in cel file
normal_cels='/restricted/projectnb/camplab/home/jia_liu/stemaway_folder/stemaway/unzipdata'
affy_norm = ReadAffy(celfile.path=normal_cels)

#normalization
norm_arry <- rma(affy_norm)
norm_matrix=exprs(norm_arry)
matrix=norm_matrix[,1:98]

#read in metadata
meta=as.data.frame(metadata)
View(meta)

#batch correction 
library(sva)
model <- model.matrix(~1, data = meta)
combat_norm_matrix <- ComBat(dat = matrix, batch = meta$BATCH, mod = model)

##visualization to compare data before and after preprocessing

#pca before normalization  
library(ggplot2)
library(grid)
library(gridExtra)

before=exprs(affy_norm)
before=before[,1:98]
#matrix
#before=matrix
colnames(before) <- factor(meta$CN)
ta=t(before)
ta=log(ta)
df_pca <- prcomp(ta,center = F, scale. = F)
df_out <- as.data.frame(df_pca$x)
group <- factor(meta$CN)
sample_lab=meta$`!Sample_geo_accession`

ggplot(df_out,aes(x=PC1,y=PC2,label= sample_lab),size=5) + 
  geom_point(aes(color=group)) +geom_text(hjust=-0.1,vjust=0.1,size=2) +
  ggtitle("PCA before normalization")

#pca before batch correction (after normalization)
before=matrix
colnames(before) <- factor(meta$CN)
ta=t(before)
ta=log(ta)
df_pca <- prcomp(ta,center = F, scale. = F)
df_out <- as.data.frame(df_pca$x)
group <- factor(meta$CN)
sample_lab=meta$`!Sample_geo_accession`

ggplot(df_out,aes(x=PC1,y=PC2,label= sample_lab),size=5) + 
  geom_point(aes(color=group)) +geom_text(hjust=-0.1,vjust=0.1,size=2) +
  ggtitle("PCA before batch correction")

#pca after normalization and batch correction
after=combat_norm_matrix
colnames(after) <- factor(meta$CN)
ta=t(after)
ta=log(ta)
df_pca <- prcomp(ta,center = F, scale. = F)
df_out <- as.data.frame(df_pca$x)
group <- factor(meta$CN)
sample_lab=meta$`!Sample_geo_accession`
namelist=colnames(before)

ggplot(df_out,aes(x=PC1,y=PC2,label= sample_lab),size=1) + 
  geom_point(aes(color=group)) +geom_text(hjust=-0.1,vjust=0.1,size=2) + 
  ggtitle("PCA after batch correction") #+ geom_text_repel()

#heatmap before batch correction
library(pheatmap)
before=matrix
colnames(before) <- factor(meta$CN)
coln=colnames(before)
dismat <- 1-cor(before)
nm=colnames(before)
annotation_col = data.frame(SampleType = nm)
rownames(annotation_col) = sample_lab
colnames(dismat) = sample_lab
pheatmap(dismat, annotation_col = annotation_col,show_rownames=F,main = 'hierarchical clustering heatmap before batch correction (after normalization)')

#heatmap after
coln=colnames(after)
dismat <- 1-cor(after)
nm=colnames(after)
annotation_col = data.frame(SampleType = nm)
rownames(annotation_col) = sample_lab
colnames(dismat) = sample_lab
pheatmap(dismat, annotation_col = annotation_col,show_rownames=F,main = 'hierarchical clustering heatmap after batch correction')




