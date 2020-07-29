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



#DE analysie
##################ANNOTATION
library('hgu133plus2.db')
x <- hgu133plus2.db
keys(x)
columns(x)
keytypes(x)

# Map Probe IDs to Gene Symbols
cols <- c("SYMBOL", "GENENAME")
Sym_ID2 <- AnnotationDbi::select(x, keys=as.character(rownames(combat_norm_matrix)), columns='SYMBOL')
table=na.omit(Sym_ID2)

# make a dataframe contains gene expression, probeid and annotation
ind <- match(rownames(combat_norm_matrix), table$PROBEID)
probeid=rownames(combat_norm_matrix)
ge <- as.data.frame(combat_norm_matrix)
ge$probeid=probeid
ge$Gene_Symbol <- table$SYMBOL[ind]

# drop duplicated rows
annot_gene<-ge[!duplicated(ge$Gene_Symbol), ] 
#annot_gene<-ge[unique(ge$Gene_Symbol), ] #you can also use unique() to drop duplication


#######filtering lowest 4% genes
aa=annot_gene[,1:98]
annot_gene$gesum=rowSums(aa)
#library(genefilter)
lowthre=quantile(annot_gene$gesum,0.04)
filt_gene=annot_gene[which(annot_gene$gesum>lowthre),]
filt_gene=na.omit(filt_gene)
filt_matrix=filt_gene[,1:98]
rownames(filt_matrix)<-filt_gene$Gene_Symbol
#rownames(filt_matrix)<-rownames(filt_gene)
#######use limma to do filtering
library(limma)
group <- factor(meta$CN,levels = c("Normal","Cancer"),ordered = F)
design <- model.matrix(~group)
colnames(design) <- levels(mygroup)
rownames(design) <- colnames(filt_matrix)
fit <- lmFit(filt_matrix, design)
fit <- eBayes(fit)
output <- topTable(fit, coef=2,n=Inf)
#sum(output$adj.P.Val<0.05)
sigout=output[output$adj.P.Val<0.05,]
# I specify coef=2 because we are interested in the difference between groups, not the intercept.

############################## volcano plot
library(EnhancedVolcano)
## Sort by ordered adj.P.Val
res_tableOE_ordered <- output[order(abs(output$logFC),decreasing = TRUE), ] 
## Create a column to indicate which genes to label
res_tableOE_ordered$genelabels <- ""
res_tableOE_ordered$genelabels[1:50] <- rownames(res_tableOE_ordered)[1:50]
#View(res_tableOE_ordered)
EnhancedVolcano(res_tableOE_ordered,
                'logFC', y = 'adj.P.Val',
                lab = res_tableOE_ordered$genelabels)

#heatmap for significant genes
hmdata=sigout[abs(sigout$logFC)>2.8,]
library(pheatmap)
ind <- match(rownames(hmdata), rownames(filt_matrix))
myre= filt_matrix[ind,]
nm=meta$CN
annotation_col = data.frame(SampleType = nm)
rownames(annotation_col) = sample_lab
colnames(myre) = sample_lab
pheatmap(myre, annotation_col = annotation_col,main = 'top 50 significant DEGs')

#another way to do heatmap
top.genes <- topTable(fit, coef=2, p.value = 0.05, number=50)
heatmap_InputMatrix<-filt_matrix[(rownames(filt_matrix) %in%  rownames(top.genes)),]
colnames(heatmap_InputMatrix)=sample_lab
pheatmap(heatmap_InputMatrix, annotation_col = annotation_col,main = 'top 50 significant DEGs')


