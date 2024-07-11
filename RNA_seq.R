#for RNA seq ana

#1. DESeq2 analysis

library(DESeq2)


setwd("")
dat <- read.delim('gene_count.xls', row.names = 1, sep = '\t', check.names = FALSE)
dat <- dat[,1:6] #'gene_count.xls' contains all gene counts of 6 samples
dat <- dat[rowMeans(dat)>1,]

coldata <- data.frame(condition = factor(rep(c('20X', '1000X'), each = 3), levels = c('20X', '1000X')))
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', '20X', '1000X'))
res
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'DESeq2.xls', col.names = NA, sep = '\t', quote = FALSE)

res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),] 
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]  
res1_total <- rbind(res1_up,res1_down)

#normalization by counts {DESeq2}

normalized_counts <- counts(dds1, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing = T),]
write.table(normalized_counts, 'DESeq2_normalized.xls', row.names = T, col.names = T, sep = '\t', quote = FALSE)

#merge DEGs and normalized counts


res_nor_all <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by ="row.names", sort=FALSE)
res_nor_DEGs <- merge(as.data.frame(res1_total), as.data.frame(counts(dds1, normalized=TRUE)), by ="row.names", sort=FALSE)
res_nor_up <- merge(as.data.frame(res1_up), as.data.frame(counts(dds1, normalized=TRUE)), by ="row.names", sort=FALSE)
res_nor_down <- merge(as.data.frame(res1_down), as.data.frame(counts(dds1, normalized=TRUE)), by ="row.names", sort=FALSE)


#2. Enrichment

#Build local database by AnnotationForge and perform enrichment analysis by clusterProfiler

library(AnnotationForge) 
library(clusterProfiler) 
library(dplyr)
library(job)
library(XML)


rm(list=ls())
setwd("")

#Import the annotations file
emapper <- rio::import('Emapper.annotations.txt') #'Emapper.annotations.txt' was a reference genome annotation file annotated by http://eggnog-mapper.embl.de/
emapper[emapper==""]<-NA


#Extract GO term
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit() 
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()

#Build an empty dataframe

gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())

#Fill the gene2go dataframe

job::job({for (row in 1:nrow(gos)) { 
  the_gid <- gos[row, "query"][[1]] 
  the_gos <- str_split(gos[row,"GOs"], ",", simplify = FALSE)[[1]] 
  df_temp <- data_frame(GID = rep(the_gid, length(the_gos)), 
                        GO = the_gos, 
                        EVIDENCE = rep("IEA", length(the_gos))) 
  gene2go <- rbind(gene2go, df_temp)}})

gene2go$GO[gene2go$GO=="-"]<-NA 
gene2go<-na.omit(gene2go)


#Information sets of the reference genome 
tax_id = "1196095 "
genus = "Gilliamella"
species = "apicola"

#duplicated
gene2go <- unique(gene2go) 
gene2go <- gene2go[!duplicated(gene2go),] 
gene2ko <- gene2ko[!duplicated(gene2ko),] 
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]

#Build an OrgDb database
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               maintainer='CF Yang <cfyang07@126.com>',
               author='CF Yang <cfyang07@126.com>',
               pathway=gene2pathway,
               version="1.32.0",
               outputDir = "..",
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")

ricenew.orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")

#Import OrgDb database

install.packages("../org.Gapicola.eg.db", repos=NULL, type="sources")
library(org.Gapicola.eg.db)
columns(org.Gapicola.eg.db)
keys(org.Gapicola.eg.db)


#GO enrichment

ego <- enrichGO(gene = res_nor_DEGs,
                OrgDb = org.Gapicola.eg.db,
                keyType = "GID",
                ont = "all",
                qvalueCutoff = 0.05,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")

write.table(ego, file = "GO.xls",sep = "\t", quote = F)
view(ego)

