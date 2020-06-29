#week4 DE analysie
##################ANNOTATION
library('hgu133plus2.db')
x <- hgu133plus2.db
keys(x)
columns(x)
keytypes(x)

# Map Probe IDs to Gene Symbols
cols <- c("SYMBOL", "GENENAME")
Sym_ID2 <- AnnotationDbi::select(x, keys=as.character(rownames(affy_matrix)), columns='SYMBOL')
table=na.omit(Sym_ID2)

# make a dataframe contains gene expression, probeid and annotation
ind <- match(rownames(affy_matrix), table$PROBEID)
probeid=rownames(affy_matrix)
ge <- as.data.frame(affy_matrix)
ge$probeid=probeid
ge$Gene_Symbol <- table$SYMBOL[ind]

# drop duplicated rows
annot_gene<-ge[!duplicated(ge$Gene_Symbol), ] 
#annot_gene<-ge[unique(ge$Gene_Symbol), ] #you can also use unique() to drop duplication


#######filtering lowest 4% genes
aa=annot_gene[,1:34]
annot_gene$gesum=rowSums(aa)
library(genefilter)
lowthre=quantile(annot_gene$gesum,0.04)
filt_gene=annot_gene[which(annot_gene$gesum>lowthre),]
filt_gene=na.omit(filt_gene)
filt_matrix=filt_gene[,1:34]
rownames(filt_matrix)<-filt_gene$Gene_Symbol

#######use limma to do filtering
library(limma)
group <- factor(c(rep("cancer",17), rep("control",18)))
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design) <- colnames(filt_matrix)
fit <- lmFit(filt_matrix, design)
fit <- eBayes(fit)
output <- topTable(fit, coef=2,n=Inf)
#sum(output$adj.P.Val<0.05)
sigout=output[output$adj.P.Val<0.05,]

# I specify coef=2 because we are interested in the difference between groups, not the intercept.
t <- topTable(fit, coef=2, n=nrow(my.second.sub), adjust='BH') #the report, adjustment is optional

############################## volcano plot
library(ggplot2)
threshold_OE <- sigout$adj.P.Val < 0.05
length(which(threshold_OE))
## Add logical vector as a column (threshold) to the res_tableOE
sigout$threshold <- threshold_OE 
##set gene label
## Sort by ordered padj
res_tableOE_ordered <- sigout[order(sigout$adj.P.Val), ] 
## Create a column to indicate which genes to label
res_tableOE_ordered$genelabels <- ""
res_tableOE_ordered$genelabels[1:10] <- rownames(res_tableOE_ordered)[1:10]
View(res_tableOE_ordered)
## Volcano plot
ggplot(res_tableOE_ordered) +
  geom_point(aes(x=logFC, y=-log10(adj.P.Val), color = ifelse(abs(logFC)>1.5&abs(-log10(adj.P.Val))>5,"red","grey"))) +
  ggtitle("DE volcano plot") +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = genelabels)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_vline(
    xintercept = c(-1.5,1.5),
    col = "black",
    linetype = "dotted",
    size = 1) +
  geom_hline(
    yintercept = 5,
    col = "black",
    linetype = "dotted",
    size = 1)+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



library(EnhancedVolcano)
##set gene label
## Sort by ordered adj.P.Val
res_tableOE_ordered <- sigout[order(sigout$adj.P.Val), ] 
## Create a column to indicate which genes to label
res_tableOE_ordered$genelabels <- ""
res_tableOE_ordered$genelabels[1:10] <- rownames(res_tableOE_ordered)[1:10]
View(res_tableOE_ordered)
s_genes <- c('abo', 'Ace')
EnhancedVolcano(res_tableOE_ordered,
                'logFC', y = 'adj.P.Val',
                lab = res_tableOE_ordered$genelabels)
