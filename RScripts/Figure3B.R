#rm(list=ls())
#setwd("~/Documents/my_tools/fungin/1")

###########Figure 3B##################
base::source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helper_functions.fungin_new.R")
## Get some mutation data
library(maftools)
####MAF Object with both SNVs and CNVs##############
o<-load(file = paste0("/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/driverGenesByCancer.rda"))
#CNVs<-rbind(l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,l1$`Neuroendocrine-lung`)
CNVs<-l1$Pheochromocytoma
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
custom.cn.data<-custom.cn.data[!(custom.cn.data$Gene %in% c("HLA-A","HLA-B","HLA-C")),]
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
wesSNVs<-l1$Pheochromocytoma
wesSNVs@data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#wesSNVs@clinical.data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@clinical.data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#head(custom.cn.data)
maf.plus.cn = read.maf(maf = wesSNVs@data,cnTable = custom.cn.data,verbose = FALSE)
mafObject<-maf.plus.cn
##################################################################################

mutatedGenes<-mafObject@data$Hugo_Symbol
RTNo<-as.character(mafObject@data$Tumor_Sample_Barcode)
mutRTMap<-data.frame(Gene=mutatedGenes,RTNo=RTNo)

o1<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
tpmValuesFull<-tpmValuesFull1
tpmValuesFull$GeneName<-row.names(tpmValuesFull)#mapping$GeneName
cancer <- "Pheochromocytoma"#"Neuroendocrine"#"Medullary thyroid carcinoma"#"Adrenocortical carcinoma"#
###Get the cancer-specific genes from the RData
cancerSpecificGenes1<-(cellSpecificGenesOutput %>% dplyr::filter(Group %in% "Tissue-Enriched") %>% dplyr::filter(grepl(cancer,Tissue)))$Gene
###Genes from co-expressed gene modules
#cancerSpecificGenes1<-c(cancerSpecificGenes1,scan("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/datKME/module_MM.red.txt",character()))
#cancerSpecificGenes1<-mapping[cancerSpecificGenes,]$GeneName
#cancerSpecificGenes1<-c()
##############################################
#finalTPM<- tpmValues %>% dplyr::filter(GeneName %in% mutatedGenes)
finalTPM<- tpmValuesFull %>% dplyr::filter(GeneName %in% c(mutatedGenes,cancerSpecificGenes1))
expData <- data.frame(do.call("rbind",lapply(unique(finalTPM$GeneName),FUN = function(col) colMeans(finalTPM[finalTPM$GeneName == col,c(seq(1, (ncol(finalTPM) - 1)))]))))
row.names(expData)<-unique(finalTPM$GeneName)
#####Calculate the Z-Score for each gene##
r<-data.frame(t(apply(expData, 1, scale)))
colnames(r)<-colnames(expData)
expData<-r
######################################
colnames(expData) <- gsub("^X","",colnames(expData))
expDataCancer<-expData[,row.names(features[grepl(cancer,features$New.Diagnosis),])]
expDataRTNo<-unlist(lapply(colnames(expDataCancer), function(x){unlist(str_split(x,pattern = "_"))[2]}))
colnames(expDataCancer)<-expDataRTNo
##Get the expression from the mutated sample
expSample<-unlist(lapply(row.names(expDataCancer), function(x){
  RTNo<-unlist(mutRTMap %>% dplyr::filter(Gene %in% x) %>% dplyr::select(RTNo))
  c<-unlist(expDataCancer[x,(colnames(expDataCancer) %in% RTNo)])
  if(length(c) > 0)
  {
    return(mean(c,na.rm=T))
  }else
  {
    return(NA)
  }
}))
names(expSample)<-row.names(expDataCancer)
##Get the expression from the all the samples for genes
# expSample<-unlist(lapply(row.names(expDataCancer), function(x){
#   #RTNo<-unlist(mutRTMap %>% dplyr::filter(Gene %in% x) %>% dplyr::select(RTNo))
#   c<-unlist(expDataCancer[x,])
#   if(length(c) > 0)
#   {
#     return(median(c,na.rm=T))
#   }else
#   {
#     return(NA)
#   }
# }))
# names(expSample)<-row.names(expDataCancer)
##Get the expression from the cancer-specific gene sample
expRNA<-unlist(lapply(names(expSample[is.na(expSample)]), function(x){
  c<-unlist(expDataCancer[x,])
  return(median(c,na.rm=T))
}))
names(expRNA)<-names(expSample[is.na(expSample)])
#expRNA<-c()
#########################################
#expSample<-log2(expSample + 1)
exp<-c(expSample[!is.na(expSample)],expRNA)
#exp<-log2(unlist(apply(expDataCancer,1,FUN=function(x) median(x,na.rm=T))) + 1)
genesWithNoExpData<-setdiff(c(mutatedGenes,cancerSpecificGenes1),names(exp))
##Check for TF in the data
TFList<-read.table("~/MyPART/AllSamplesPipeliner/Homo_sapiens_TF.txt",header = T,sep = "\t")
isTFVector<-c(names(exp),genesWithNoExpData) %in% TFList$Symbol
rnaseq_results <- data.frame(gene_symbol=c(names(exp),genesWithNoExpData),logFoldChange=rep(NA,length(exp) + length(genesWithNoExpData)),pvalue=rep(NA,length(exp) + length(genesWithNoExpData)),expression=c(exp,rep(NA,length(genesWithNoExpData))),isTF=isTFVector,stringsAsFactors = F)

rnaseq_results[rnaseq_results$gene_symbol %in% wesSNVs@data[]$Hugo_Symbol,]$pvalue<-"SNVs"
rnaseq_results[rnaseq_results$gene_symbol %in% CNVs[CNVs$driver %in% c("DEL"),]$gene,]$pvalue<-"CNVs-DEL"
rnaseq_results[rnaseq_results$gene_symbol %in% CNVs[CNVs$driver %in% c("AMP","PARTIAL_AMP"),]$gene,]$pvalue<-"CNVs-AMP"
rnaseq_results[rnaseq_results$gene_symbol %in% intersect(wesSNVs@data$Hugo_Symbol,CNVs$gene),]$pvalue<-"CNVs-SNVs"

###Adding annotation about the drug bank data
drugTargets<-read.table(paste0("/Users/jaina13/myPART/AllSamplesPipeliner/","CancerSpecificGenes/DrugBankTargets-Format.tsv"),header = TRUE,sep = "\t")
drugBankGeneMapping<-read.table(paste0("/Users/jaina13/myPART/AllSamplesPipeliner/","CancerSpecificGenes/genesDrugBank.tsv"),header = TRUE,sep = "\t")
drugTargetsMap<-merge(x=drugTargets,y=drugBankGeneMapping,by.x="targetID",by.y="gene_claim_name",all.x=FALSE,all.y=FALSE)
rnaseq_results$isdrugTarget<-"N"
rnaseq_results$isdrugTarget[rnaseq_results$gene_symbol %in% unique(drugTargetsMap$gene_name)]<-"Y"

#### There are basically two ways to use this framework:
# 1. Start with a large list of genes (say 50-200?) and make a network out of the connections among them, or
# 2. Start with a small list of genes (say < 10), and grow a network out of their interacting neighbors

## Large list of genes by taking union of top genes from expression and mutation data
#query_genes_large <- union(laml@gene.summary$Hugo_Symbol[1:top_n_mutated_genes], dummy_rnaseq_results$gene_symbol)
query_genes_large <- c(mutatedGenes,cancerSpecificGenes1)#union(laml@gene.summary$Hugo_Symbol[1:top_n_mutated_genes], dummy_rnaseq_results$gene_symbol)

## Small list of genes
#query_genes_small <- laml@gene.summary$Hugo_Symbol[1:5]

# hubGenes<-analyzePPINetwork(unique(query_genes_large),paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WES-CNVs-CGGenes-Pheochromocytoma.STRING.PPI.png"),ppiEdgeThreshold = 700)
# write.table(hubGenes,paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WES-CNVs-CGGenes-Pheochromocytoma.hubGenes.PPI.txt"),quote = F,row.names = F,sep = "\t")


#########################
#######  StringDB #######
#########################
network_source="stringDB"
network_data_dir=file.path("/Users/jaina13/myPART/WGSCode/data")
outdir=file.path("/Users/jaina13/myPART/WGSCode/data/figures", network_source)
#if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

## StringDB's network edges have a "combined_score" attribute that ranges from 0-999
#string_score_threshold=400
string_score_threshold=700
fullgraph <- fungin_initialize(fungin_data_dir = network_data_dir, sourcedb = network_source, genome = "hg38", string_score_threshold = string_score_threshold)
fungin_query <- fungin_trim(fullgraph, query_genes = query_genes_large, get_neighbors = F,get_largest_connected_subnet = F)

# 1. Make plot with MAF data only
# fungin_anno <- fungin_annotate(fungin_query, maf=mafObject, min_mutated_samples = 1)
#
# mynetworkplot <- fungin_plot(fungin_anno, fillVar = "Nmutated")
# networkplot_image <- file.path(outdir, "ex1.maf_only.pdf")
# ggsave(filename = networkplot_image, plot=mynetworkplot, width=8, height=6)

# 2. Make plot with MAF and RNA data
## The 'diff_exp_results' argument is a data frame with at lest three columns: gene, logFC, and pval
## These column names can be set with the *_column arguments (the defaults are set for topTable() output from limma)
base::source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helper_functions.fungin_new.R")
fungin_anno <- fungin_annotate(fungin_query, maf=mafObject, min_mutated_samples = 1,
                               diff_exp_results = rnaseq_results,
                               gene_column = "gene_symbol", fc_column = "logFoldChange", pval_column = "pvalue",expression_column="expression")

mynetworkplot <- fungin_plot(fungin_anno,fillVar = "Expression")
#mynetworkplot <- mynetworkplot + theme(plot.margin = unit(c(1,1,1,1), "cm")) +
#mynetworkplot_int <-fungin_plot_interactive(fungin_anno,fill_var = "Expression")
ggsave(filename = paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/Pheochromocytoma_PPI_ZScore_csGenes_700.pdf"), plot=mynetworkplot, width=12, height=10)
# networkplot_image <- file.path(outdir, "maf_with_rna_ZScore.pdf")
# ggsave(filename = paste0(mywd,"PPINetwork/ACC_PPI_ZScore_csGenes_900.pdf"), plot=mynetworkplot, width=12, height=8)

####Reduce the density of the network using the igraph function
test.layout <- layout_(fungin_anno,with_dh(weight.edge.lengths = edge_density(fungin_anno)/1000))
plot(fungin_anno, layout = test.layout)

# 3. Make plot with RNA data only, and growing the network by including neighbors
# fungin_query <- fungin_trim(fullgraph, query_genes = query_genes_small, get_neighbors = T)
# fungin_anno <- fungin_annotate(fungin_query, maf=laml, min_mutated_samples = 3,
#                                diff_exp_results = dummy_rnaseq_results,
#                                gene_column = "gene_symbol", fc_column = "logFoldChange", pval_column = "pvalue")
#
# mynetworkplot <- fungin_plot(fungin_anno)
# networkplot_image <- file.path(outdir, "ex1.rna_only.pdf")
# ggsave(filename = networkplot_image, plot=mynetworkplot, width=8, height=6)
