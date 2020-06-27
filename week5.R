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

#####################2:富集分析的背景基因集
gene.df <- bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "GO"), OrgDb = org.Hs.eg.db)
head(gene.df,2)

######################3.1: GO analysis
ggo <- groupGO(gene = gene.df$ENTREZID,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)

#####################3.2 enrichGo analysis
ego_ALL <- enrichGO(gene = gene.df$ENTREZID, 
                    #背景基因集
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) #Gene ID 转成gene Symbol ，易读
head(ego_ALL,2)

#####################3.2.1 setReadable函数进行转化
ego_MF <- enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
#####################3.3 GSEA分析 （暂略）

#####################3.4 write out
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)
##可视化--点图
dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
##可视化--条形图
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term
##可视化--
plotGOgraph(ego_MF)

#####################4: KEGG分析
#####################4.1 候选基因进行通路分析
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
####################write out
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)
dotplot(kk,title="Enrichment KEGG_dot")


########################5: 注释文件、注释库
#####5.1 待富集的基因list
deg <- names(geneList)[abs(geneList)>2] 


########### WikiPathways analysis
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

aa=sigout$logFC
bb=as.character(gene.df$ENTREZID)
names(aa) <- bb
aa <- sort(aa, decreasing = TRUE)
ewp2 <- GSEA(aa, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)


