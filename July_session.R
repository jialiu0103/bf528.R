#read in cel file
normal_cels='/stemaway/unzipdata'
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




####################################### functional analysis ##################################################
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
keytypes(org.Hs.eg.db) 
####################1:prepare dataset
## assume that 1st column is ID, 2nd column is fold change
## feature 1: numeric vector
geneList <- sigout$logFC
#geneList <- data[,2]

## feature 2: named vector
aa=rownames(sigout)
bb=as.character(aa)
names(geneList) <- bb

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

##Suppose we define fold change greater than 2 as DEGs:
gene <- names(geneList)[abs(geneList) > 2]

#####################2:convert gene symbol to geneid
gene.df <- bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "GO"), OrgDb = org.Hs.eg.db)
head(gene.df,2)

######################3.1: GO analysis
ggo <- groupGO(gene = gene.df$ENTREZID,keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
aa=setReadable(ggo, OrgDb = org.Hs.eg.db)
barplot(aa)
#####################3.2 enrichGo analysis
ego_ALL <- enrichGO(gene = bb, 
                    OrgDb = org.Hs.eg.db, #human: OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #or CC,  BP,  MF
                    pAdjustMethod = "BH", #could be "holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
                    pvalueCutoff = 1, #P value filtering
                    qvalueCutoff = 1,
                    readable = TRUE) #Gene ID convert to gene Symbol, make it readable 
head(ego_ALL,2)

#####################3.2.1 setReadable to make it readable 
ego_MF <- enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
#####################3.3 GSEA

#####################3.4 write out
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)
#dotplot
dotplot(ego_MF,title="EnrichmentGO_MF_dot")#dotplot
##barplot
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")#barplot, first 20 Term
##graph
plotGOgraph(ego_MF)

#####################4: KEGG analysis
#####################4.1 pathway analysis 
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = 'hsa', #KEGG can use organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
####################write out
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)
dotplot(kk,title="Enrichment KEGG_dot")


########################5: annotate 
#####5.1 set a gene list
deg <- names(geneList)[abs(geneList)>2] 


####################### WikiPathways analysis
##WikiPathways is a continuously updated pathway database curated by a community of researchers and pathway 
##enthusiasts. WikiPathways produces monthly releases of gmt files for supported organisms at data.wikipathways.org.
library(magrittr)
library(tidyr)

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wpgene <- read.gmt(wpgmtfile)
wp2gene <- wpgene

path='/restricted/projectnb/camplab/home/jia_liu/stemaway_folder/stemaway/GSE21510_RAW/geneset.gmt'
wpgmtfile <- read.gmt(path)
wp2gene <- wpgmtfile

wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)



aa=setReadable(ewp, OrgDb = org.Hs.eg.db, keyType='ENTREZID')
barplot(aa)


aa=sigout$logFC
bb=as.character(gene.df$ENTREZID)
names(aa) <- bb
aa <- sort(aa, decreasing = TRUE)
ewp2 <- GSEA(aa, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
bb=setReadable(ewp2, OrgDb = org.Hs.eg.db, keyType='ENTREZID')
saveRDS(bb)


###############cnetplot
myde<- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb = org.Hs.eg.db)
mydee<-myde$ENTREZID
myedo <- enrichDGN(mydee)
myedox <- setReadable(myedo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)
cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
