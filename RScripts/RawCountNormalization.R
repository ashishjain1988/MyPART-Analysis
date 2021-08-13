library(DESeq2)
library(apeglm)
library(edgeR)
library(WebGestaltR)

source("~/myPART/MyPART-Analysis/RScripts/helperFunctions.R")

##Path to the base directory with the data
parent<-"~/myPART/AllSamplesPipeliner/"
countsFile <- "RawCountFile_RSEM_genes.txt"
counts <- read.table(paste0(parent,countsFile),header=TRUE, row.names = 1)
colnames(counts) <- gsub("^Sample_","",colnames(counts))
colnames(counts) <- unlist(lapply(colnames(counts),function(x){i<-unlist(gregexpr(pattern ="_S",x));return(substr(x,0,i[length(i)]-1));}))

##Features to make design
features<-read.table(paste0("~/myPART/features.txt"),header = TRUE,sep = "\t",row.names = 1)
counts<-counts[,row.names(features)]

##Make subtypes of NET based on primary tissue
NETSamples<-features[features$New.Diagnosis == "Neuroendocrine",]
features[features$New.Diagnosis == "Neuroendocrine",]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)
###Removing the samples with bad quality (low TIN) from the analysis
features<-features[!(row.names(features) %in% c("47_RT00085_X007R_resent","76_RT00105_X336R")),]
features<-features[features$Tissue != "Human Universal Reference Total RNA",]

##Groups to use for DEGs Analysis
samplesToCompare<-row.names(features[features$Subtype %in% c("Endocrine") | features$New.Diagnosis %in% c("Normal"),])
#samplesToCompare<-row.names(features[features$New.Diagnosis %in% c("Neuroendocrine"),])#"Adrenocortical carcinoma",c("Gastrointestinal stromal tumor"),"Carcinoid tumor"#
# normalSamples<-row.names(features[features$Diagnosis %in% c("Normal") & features$Tissue %in% features[samplesToCompare,]$Tissue,])
# samplesToCompare<-row.names(features[features$Diagnosis %in% c("Gastrointestinal stromal tumor") & features$Tissue %in% features[normalSamples,]$Tissue,])
# samplesToCompare<-c(samplesToCompare,normalSamples)
countTable <- counts[,samplesToCompare]
featuresFilt<-features[samplesToCompare,]

#Filtering based on the CPM values
CPM_THRESHOLD <- 0.5
keep.rows <- which(rowMeans(edgeR::cpm(countTable)) >= CPM_THRESHOLD)
allcounts <- countTable[keep.rows,]

#Make Diagnosis based for different comparisons
#NETSamples<-featuresFilt[featuresFilt$PrimaryCancerTissue == "small intestine",]
#featuresFilt[featuresFilt$PrimaryCancerTissue == "small intestine",]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)
##The count from rsem is not integer so have to round it
#featuresFilt$New.Diagnosis<-as.factor(featuresFilt$New.Diagnosis)
dds <- DESeqDataSetFromMatrix(countData = round(allcounts),colData = featuresFilt,design = ~ New.Diagnosis + Tissue)
# dds <- DESeqDataSetFromMatrix(countData = round(allcounts),colData = featuresFilt,design = ~ cluster)
##Filter from variancePartition package
#isExpressed <- rowSums(fpm(dds)>1) >= 0.5 * ncol(dds)
#dds$Gender<-relevel(dds$Gender,ref = "Female")
#dds$Diagnosis<-relevel(dds$Diagnosis,ref = "Normal")
#dds$cluster<-relevel(dds$cluster,ref = paste0("Other",c))
#dds$Diagnosis<-relevel(dds$New.Diagnosis,ref = "Normal")

##Wald test
#dds <- DESeq(dds,minReplicatesForReplace=Inf)
#dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE,nsub = nrow(dds))
vsdNormalizedCounts<-assay(vsd)
save(featuresFilt,vsdNormalizedCounts,file = paste0(parent,"EndocrineSubgroupResults/vsdNormalizedCounts.rda"))
