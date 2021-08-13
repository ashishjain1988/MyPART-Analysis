library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)

#### Inputs
mywd="/Users/jaina13/myPART/FusionAnalysis/"  ## Working directory
fusion_out_dir <- "/Users/jaina13/myPART/FusionAnalysis/arriba_results/"     ## Location of annotsv output files

fusion_suffix <- "*.fusions.tsv"
sample_files <- list.files(path = fusion_out_dir,pattern = fusion_suffix,recursive = T, full.names = T)
names(sample_files) <- path_file(sample_files)
all_fusion_results <- lapply(sample_files, function(currfile) {
  samplename <-  unlist(str_split(path_file(currfile),pattern = "_S"))[1]
  # print(samplename)
  #curr_results <- read.table(currfile, header=T, sep="\t", comment.char = "", stringsAsFactors = F)
  curr_results <- fread(currfile,  sep="\t")
  if (nrow(curr_results) > 1) {
    curr_results$sample <- samplename
  } else {
    curr_results <- NA
  }
  return(curr_results)
})
all_fusion_results <- all_fusion_results[!is.na(all_fusion_results)]
sNames<-unlist(lapply(names(all_fusion_results),FUN=function(x){
  return(unlist(str_split(x,pattern = "_S"))[1])
}))
names(all_fusion_results)<-sNames

#Feature files
features.file <- "/Users/jaina13/myPART/features.txt"
features <- read.table(features.file, sep="\t", header = T, stringsAsFactors = F)
##Make subtypes of NET based on primary tissue
NETSamples<-features[features$New.Diagnosis == "Neuroendocrine",]
features[features$New.Diagnosis == "Neuroendocrine",]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)
row.names(features)<-paste0("Sample_",features$SampleID)
features<-features[names(all_fusion_results),]
###Removing the samples with bad quality (low TIN) from the analysis
features<-features[!(row.names(features) %in% c("47_RT00085_X007R_resent","76_RT00105_X336R")),]

normal_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c("Normal") ))#& Tissue %in% c("Adrenal")))
normal_fus_res <-all_fusion_results[normal_samples]
normFusionsDF <- do.call("rbind",normal_fus_res) %>% dplyr::filter(confidence == "high")
normFusionsDF$Fusion<-paste0(normFusionsDF$`#gene1`,"-",normFusionsDF$gene2,"-",normFusionsDF$breakpoint1,"-",normFusionsDF$breakpoint2)
mywd<-"/Users/jaina13/myPART/FusionAnalysis/arriba_results//"
names<-unique(features$New.Diagnosis)
l<-list()

for(name in names[2:length(names)])
{
  tumor_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c(name)))#,"Adrenocortical carcinoma""Gastrointestinal stromal tumor"
  tumor_fus_res <- all_fusion_results[tumor_samples]
  allFusionsDF <- do.call("rbind",tumor_fus_res)

  #Filtering of the fusion genes based on confidence
  allFusionsDF$Fusion<-paste0(allFusionsDF$`#gene1`,"-",allFusionsDF$gene2,"-",allFusionsDF$breakpoint1,"-",allFusionsDF$breakpoint2)
  allFusionsDFFilt<-allFusionsDF %>% dplyr::filter(confidence == "high" & (!Fusion %in% normFusionsDF$Fusion))
  #& (!X.FusionName %in% normFusionsDF$X.FusionName))
  if(nrow(allFusionsDFFilt) > 0)
  {
    proteinCodingGenes<-read.table(paste0("/Users/jaina13/myPART/FusionAnalysis/hg38.proteinCodingGenes.Ensembl.txt"),sep = "\t",header = T)
    allFusionsDFFilt_pc <- allFusionsDFFilt %>% dplyr::filter(`#gene1` %in% proteinCodingGenes$Gene.name & gene2 %in% proteinCodingGenes$Gene.name)
    if(nrow(allFusionsDFFilt_pc)>0)
    {
      uniqueFusions <- allFusionsDFFilt_pc %>% dplyr::group_by_at(vars(Fusion)) %>% dplyr::summarize(split_reads1=paste(split_reads1,collapse = ","),
                                                                                                     split_reads2=paste(split_reads2,collapse = ","),
                                                                                                     type=paste(type,collapse = ","),confidence=paste(confidence,collapse = ","),
                                                                                                     "strand1(gene/fusion)"=paste(`strand1(gene/fusion)`,collapse = ","),
                                                                                                     "strand2(gene/fusion)"=paste(`strand2(gene/fusion)`,collapse = ","),
                                                                                                     Sample = paste(sample,collapse = ","),
                                                                                                     Number = dplyr::n()) %>% arrange(desc(Number))
      name<-str_replace_all(name,pattern = "[ /]",replacement = ".")
      l[[name]]<-uniqueFusions
    }
  }
}

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(l, paste0(mywd,"unique.fusion.protein-coding.genes.arriba.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)
