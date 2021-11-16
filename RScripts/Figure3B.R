#rm(list=ls())
#setwd("~/Documents/my_tools/fungin/1")

library(stringr)
###########Figure 3B##################
base::source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helper_functions.fungin_new.R")
## Get some mutation data
library(maftools)
####MAF Object with both SNVs and CNVs##############
o<-load(file = paste0("/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/driverGenesByCancer.rda"))
CNVs<-rbind(l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,l1$`Neuroendocrine-lung`)
#CNVs<-l1$Adrenocortical.carcinoma
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
custom.cn.data<-custom.cn.data[!(custom.cn.data$Gene %in% c("HLA-A","HLA-B","HLA-C")),]
#o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/merged_somatic_variants/maf/SNVsResults-consensous.rda"))
wesSNVs<-l1$Neuroendocrine.tumor
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
cancer <- "Neuroendocrine"#"Adrenocortical carcinoma"#"Pheochromocytoma"#"Medullary thyroid carcinoma"#
###Get the cancer-specific genes from the RData
cancerSpecificGenes1<-(cellSpecificGenesOutput %>% dplyr::filter(Group %in% "Tissue-Enriched") %>% dplyr::filter(grepl(cancer,Tissue)))$Gene
csGeneMap<-(cellSpecificGenesOutput %>% dplyr::filter(Group %in% "Tissue-Enriched") %>% dplyr::filter(grepl(cancer,Tissue)))
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
  cancerType<-csGeneMap[csGeneMap$Gene %in% x,]$Tissue
  if(length(cancerType) != 0)
  {
    c<-unlist(expDataCancer[x,features[features$New.Diagnosis %in% cancerType,]$RTNo])
  }else
  {
    c<-unlist(expDataCancer[x,])
  }
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

# source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helperFunctions-WGCNA.R")
# p<-"/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/"
# g<-getGOEnrichmentWebGestaltRWGCNA(query_genes_large,"geneontology_Biological_Process_noRedundant",stringr::str_replace_all(cancer,pattern = " ",replacement = "_"),p)
# g<- getGOEnrichmentWebGestaltRWGCNA(query_genes_large,"pathway_Reactome",stringr::str_replace_all(cancer,pattern = " ",replacement = "_"),p)

## Small list of genes
#query_genes_small <- laml@gene.summary$Hugo_Symbol[1:5]

# hubGenes<-analyzePPINetwork(unique(query_genes_large),paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WES-CNVs-CGGenes-Pheochromocytoma.STRING.PPI.png"),ppiEdgeThreshold = 700)
# write.table(hubGenes,paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WES-CNVs-CGGenes-Pheochromocytoma.hubGenes.PPI.txt"),quote = F,row.names = F,sep = "\t")

##Adding pathways to the integarted PPI Network
# library(msigdbr)
# # m_df = msigdbr(species = "Homo sapiens", category = "H")
# # m_df = m_df[grepl("KRAS|CATENIN", m_df$gs_name),]
# m_df = msigdbr(species = "Homo sapiens", category="C2",subcat="REACTOME")
# #m_df = m_df[grepl("REACTOME_SIGNALING_BY_WNT_IN_CANCER", m_df$gs_name),]
# pathwayGeneCount<-m_df %>% dplyr::group_by_at(vars(gs_name)) %>% dplyr::summarize(gs_name=unique(gs_name),NoOfGenes = dplyr::n()) %>% arrange(desc(NoOfGenes))
#
# m_df = m_df[m_df$gene_symbol %in% query_genes_large,]
# pathwayGeneCount <- pathwayGeneCount[pathwayGeneCount$gs_name %in% m_df$gs_name,] %>% dplyr::arrange(gs_name)
# topPathways<-(m_df %>% dplyr::group_by_at(vars(gs_name)) %>% dplyr::summarize(gs_name=unique(gs_name),NoOfGenes = dplyr::n()) %>% arrange(desc(gs_name)))
# foldChangePath<-topPathways$NoOfGenes/pathwayGeneCount$NoOfGenes
# names(foldChangePath)<-topPathways$gs_name
# top10<-names(sort(foldChangePath,decreasing = TRUE))[1:10]
#
# m_df = m_df[m_df$gs_name %in% top10,]
#
# genepathmap <- unique(m_df[,c("gene_symbol","gs_name")])
# genepathmap[,2] <- stringr::str_wrap(gsub("_"," ",gsub("^HALLMARK ","",genepathmap[[2]])), width = 12)

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
mynetworkplot
#mynetworkplot <- mynetworkplot + theme(plot.margin = unit(c(1,1,1,1), "cm")) +
#mynetworkplot_int <-fungin_plot_interactive(fungin_anno,fill_var = "Expression")
#ggsave(filename = paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/NET_PPI_ZScore_csGenes_700-Pathway.pdf"), plot=mynetworkplot, width=16, height=10)
ggsave(filename = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/NET_PPI_ZScore_csGenes_700.pdf"), plot=mynetworkplot, width=12, height=12)
# networkplot_image <- file.path(outdir, "maf_with_rna_ZScore.pdf")
# ggsave(filename = paste0(mywd,"PPINetwork/ACC_PPI_ZScore_csGenes_900.pdf"), plot=mynetworkplot, width=12, height=8)

#MTCIntegratedGenes<-V(fungin_anno)
#PCCIntegratedGenes<-V(fungin_anno)
write.table(c(names(MTCIntegratedGenes),names(PCCIntegratedGenes)),"~/MTC-PCC.PPIGenes.txt",quote = F,col.names = F,row.names = F)

###PCSF package
library("PCSF")
data("STRING")
#data("Tgfb_phospho")
#terminals <- Tgfb_phospho
strV115<-read.table("~/Downloads/9606.protein.links.v11.5.txt",header = T,sep = " ")
geneMap<-read.table("~/Downloads/9606.protein.info.v11.5.txt",header = F,sep = "\t")
d<-merge(strV115,geneMap,by.x = "protein1",by.y = "V1",all.x = T)[,c(1:4)]
d1<-merge(d,geneMap,by.x = "protein2",by.y = "V1",all.x = T)[,c(1:5)]
strV115PPI<-d1[,c(4:5,3)]
strV115PPIFilt <- strV115PPI[as.logical(apply(strV115PPI,1,function(x){min(!is.na(x))})),]
strV115PPIFilt<-strV115PPIFilt[strV115PPIFilt$combined_score >=700,]
##Computing cost of the edges
strV115PPIFilt$combined_score <- 1-(strV115PPIFilt$combined_score/1000)
ppi <- construct_interactome(strV115PPIFilt)

terminals<-rep(1,length(rnaseq_results$gene_symbol))#abs(rnaseq_results$expression)
names(terminals)<-rnaseq_results$gene_symbol
#terminals<-terminals[!is.na(terminals)]
subnet <- PCSF((ppi), terminals, w = 2, b = 1, mu = 0.0005)
res <- enrichment_analysis(subnet)
p<-plot.PCSFe(res$subnet)
visSave(p,paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/MTC_Int.PCSF.STRINGV11.5.html"), selfcontained = TRUE, background = "white")

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

####Make Sample specific heatmap
ppiGenes<-inputGenes<-scan("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/MTC-PCC.PPIGenes.txt",character())
base::source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helper_functions.fungin_new.R")
## Get some mutation data
library(maftools)
####MAF Object with both SNVs and CNVs##############
o<-load(file = paste0("/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/driverGenesByCancer.rda"))
CNVs<-rbind(l1$Medullary.thyroid.carcinoma,l1$Pheochromocytoma)
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
custom.cn.data<-custom.cn.data[!(custom.cn.data$Gene %in% c("HLA-A","HLA-B","HLA-C")),]
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
wesSNVs<-merge_mafs(list(l1$Medullary.thyroid.carcinoma,l1$Pheochromocytoma))
wesSNVs@data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#wesSNVs@clinical.data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@clinical.data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#head(custom.cn.data)
maf.plus.cn = read.maf(maf = wesSNVs@data,cnTable = custom.cn.data,verbose = FALSE)
mafObject<-maf.plus.cn
##################################################################################
mutatedGenes<-mafObject@data$Hugo_Symbol
RTNo<-as.character(mafObject@data$Tumor_Sample_Barcode)
mutationType <-as.character(mafObject@data$Variant_Type)
#mutRTMap<-data.frame(Gene=mutatedGenes,RTNo=RTNo)
mutationInformation<-data.frame(Genes=mutatedGenes,RTNo = as.character(mafObject@data$Tumor_Sample_Barcode),MutType=mutationType)
mutationInformationFilt<-mutationInformation[mutationInformation$Genes %in% ppiGenes,]
mutationInformationFilt<-data.frame(mutationInformationFilt %>% tidyr::pivot_wider(names_from = "RTNo", values_from = "MutType"))

####Gene Expression Information of the PPI Genes#########
o1<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
tpmValuesFull<-tpmValuesFull1
tpmValuesFull$GeneName<-row.names(tpmValuesFull)#mapping$GeneName
# cancer <- "Medullary thyroid carcinoma|Pheochromocytoma"
# ###Get the cancer-specific genes from the RData
# cancerSpecificGenes1<-(cellSpecificGenesOutput %>% dplyr::filter(Group %in% "Tissue-Enriched") %>% dplyr::filter(grepl(cancer,Tissue)))$Gene
##############################################
#finalTPM<- tpmValues %>% dplyr::filter(GeneName %in% mutatedGenes)
finalTPM<- tpmValuesFull %>% dplyr::filter(GeneName %in% ppiGenes)
expData <- data.frame(do.call("rbind",lapply(unique(finalTPM$GeneName),FUN = function(col) colMeans(finalTPM[finalTPM$GeneName == col,c(seq(1, (ncol(finalTPM) - 1)))]))))
row.names(expData)<-unique(finalTPM$GeneName)
#####Calculate the Z-Score for each gene##
r<-data.frame(t(apply(expData, 1, scale)))
colnames(r)<-colnames(expData)
expData<-r
######################################
colnames(expData) <- gsub("^X","",colnames(expData))
annotData<-features[grepl("Pheochromocytoma|Medullary thyroid carcinoma",features$New.Diagnosis),]
row.names(annotData)<-annotData$RTNo
expDataCancer<-expData[,row.names(features[grepl("Pheochromocytoma|Medullary thyroid carcinoma",features$New.Diagnosis),])]
expDataRTNo<-unlist(lapply(colnames(expDataCancer), function(x){unlist(str_split(x,pattern = "_"))[2]}))
colnames(expDataCancer)<-expDataRTNo
expDataCancer$Genes<-row.names(expDataCancer)

noGenesInSNV<-setdiff(expDataCancer$Genes,mutationInformationFilt$Genes)
more.rows <- data.frame(Genes=noGenesInSNV, stringsAsFactors=F)
mutationInformationFilt[(nrow(mutationInformationFilt) + 1):(nrow(mutationInformationFilt) + nrow(more.rows)), names(more.rows)] <- more.rows
row.names(mutationInformationFilt)<-mutationInformationFilt$Genes

noSampleInSNV <- setdiff(colnames(expDataCancer),colnames(mutationInformationFilt))
mutationInformationFilt$RT00128<-rep(NA,nrow(mutationInformationFilt))
mutationInformationFilt$RT00091<-rep(NA,nrow(mutationInformationFilt))


##Final Data Frame###
annotData<-annotData[expDataRTNo,]
annotData$Tissue<-stringr::str_to_title(annotData$Tissue)
annotData$TumorSite<-stringr::str_to_title(annotData$TumorSite)
expDataCancer<-expDataCancer[,expDataRTNo]
mutationInformationFilt<-mutationInformationFilt[row.names(expDataCancer),expDataRTNo]


source("~/myPART/MyPART-Analysis/RScripts/helperFunctions.R")
source("~/myPART/MyPART-Analysis/RScripts/Figure3A-helper.R")
hm_data<-as.matrix(expDataCancer)
column_annotations = HeatmapAnnotation(df = annotData[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = getAgeColor(annotData$Age)#colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","blue"))
col_tissue <-tissueColor[stringr::str_to_title(annotData$Tissue)]
col_disease <-diagnosisColorRNASeq[annotData$New.Diagnosis]
#tumorSiteColor <-

df<-annotData[,c("New.Diagnosis","Tissue","Age","Sex","Race","TumorSite")]
colnames(df)<-c("Diagnosis","Tissue","Age","Sex","Race","TumorSite")
ha = HeatmapAnnotation(
  df = df,
  col = list(Diagnosis = col_disease,#c("Normal" = "black","Adrenocortical carcinoma" = "red"),#col_cluster,#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             Tissue = col_tissue,
             Age = col_age,
             Sex = genderColor,#c("Male"="turquoise","Female"="brown"),
             #RIN = col_rin,
             Race =raceColor,#c("White"="#B35806","Unknown"="#FDBC6B","Other"="black","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue"),
             TumorSite=tumorSiteColor#c("Metastasis"="turquoise","Primary Site"="brown","Recurrence"="black")
             #cluster = col_cluster
  ),
  annotation_name_gp =  gpar(fontsize = 14,fontface = 2),
  #annotation_legend_param = list(title_gp = gpar(fontsize = 15, fontface = 2),labels_gp = gpar(fontsize = 12)),
  annotation_legend_param=list(labels_gp = gpar(fontsize = 12),title_gp = gpar(fontsize = 15, fontface = 2),
                               #grid_height=unit(3,"mm"),
                               #grid_width=unit(3,"mm"),
                               # ncol = 5,
                               nrow = 5,
                               legend_direction = "vertical")
)

cell_fun_custom = function(j, i, x, y, w, h, col) {
  mutationType<-as.matrix(mutationInformationFilt)[i,j]
  if(!is.na(mutationType) && mutationType != "NA")
  {
    if(mutationType == "CNV")
    {
      # grid.points(x, y,
      #             pch = 16, gp = gpar(col = "green"), size = unit(4, "mm"))
      grid.rect(x, y, w, h,
                gp = gpar(col = "green", fill = NA, lwd=5))
    }else
    {
      # grid.points(x, y,
      #             pch = 16, gp = gpar(col = "blue2"), size = unit(4, "mm"))
      grid.rect(x, y, w, h,
                gp = gpar(col = "blue2", fill = NA, lwd=5))
    }
  }
}

cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(hm_data), na.rm = T), max(abs(hm_data), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    cell_fun =  cell_fun_custom,
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=FALSE,
                                    #column_names_rot = 45,
                                    cluster_rows = TRUE,
                                    cluster_columns = TRUE,
                                    clustering_distance_columns = "euclidean",
                                    clustering_method_columns = "complete",
                                    row_names_gp = gpar(fontsize = 14,fontface = 2),
                                    #column_names_gp = gpar(fontsize = 14,fontface = 2),
                                    show_heatmap_legend = F,
                                    #heatmap_height = unit(100,"mm")
) # Turning off to control the placement

pdf(paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/MTC-PCC.PPIClustering.pdf"),width = 10,height = 10)
ComplexHeatmap::draw(cheatmap, show_annotation_legend = TRUE)
dev.off()
