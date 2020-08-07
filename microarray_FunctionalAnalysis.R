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
aa <- rownames(sigout)
bb <- as.character(aa)
names(geneList) <- bb

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

##Suppose we define fold change greater than 2 as DEGs:
gene <- names(geneList)[abs(geneList) > 2]

#####################2:convert gene symbol to geneid
gene.df <- bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "GO"), OrgDb = org.Hs.eg.db)
head(gene.df,2)

######################3.1: GO analysis
ggo <- groupGO(gene = gene.df$ENTREZID,keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = "CC",level = 5,readable = TRUE)
aa <- setReadable(ggo, OrgDb = org.Hs.eg.db)
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
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene.df$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)
aa=setReadable(ewp, OrgDb = org.Hs.eg.db, keyType='ENTREZID')
barplot(aa)


aa <- sigout$logFC
bb <- as.character(gene.df$ENTREZID)
names(aa) <- bb
aa <- sort(aa, decreasing = TRUE)
ewp2 <- GSEA(aa, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
bb <- setReadable(ewp2, OrgDb = org.Hs.eg.db, keyType='ENTREZID')
saveRDS(bb)

