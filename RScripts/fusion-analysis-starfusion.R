library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)

#### Inputs
mywd="/Users/jaina13/myPART/FusionAnalysis/"  ## Working directory
fusion_out_dir <- "/Users/jaina13/myPART/FusionAnalysis/starfusion-results"

fusion_suffix <- "*.abridged.tsv"                     ## Suffix for annotsv output files (Sample Labels = Filename - Suffix)
sample_files <- list.files(path = fusion_out_dir,pattern = fusion_suffix,recursive = T, full.names = T)
names(sample_files) <- path_file(sample_files)#basename(dirname(basename(dirname(sample_files))))
all_fusion_results <- lapply(sample_files, function(currfile) {
  # currfile=sample_files[[82]]
  samplename <-  unlist(str_split(path_file(currfile),pattern = "_S"))[1]
  # print(samplename)
  curr_results <- read.table(currfile, header=T, sep="\t", comment.char = "", stringsAsFactors = F)
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

##Fusion Gene analysis in individual groups
ffpm_cutoff=0.1
readCountCutoff<-5

normal_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c("Normal") ))#& Tissue %in% c("Adrenal")))
normal_fus_res <-all_fusion_results[normal_samples]
normFusionsDF <- do.call("rbind",normal_fus_res) %>% dplyr::filter((JunctionReadCount >= readCountCutoff | SpanningFragCount >= readCountCutoff) & FFPM >= ffpm_cutoff)
mywd<-"/Users/jaina13/myPART/FusionAnalysis/CombinedFusionGenesWONormal/"
names<-unique(features$New.Diagnosis)
l<-list()
#for(name in names)
for(name in names[2:length(names)])
{
  tumor_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c(name)))#,"Adrenocortical carcinoma""Gastrointestinal stromal tumor"
  tumor_fus_res <- all_fusion_results[tumor_samples]
  allFusionsDF <- do.call("rbind",tumor_fus_res)

  #Filtering of the fusion genes based on FFPM
  allFusionsDFFilt<-allFusionsDF %>% dplyr::filter(FFPM >= ffpm_cutoff & (JunctionReadCount >= readCountCutoff | SpanningFragCount >= readCountCutoff)
                                                   & (!X.FusionName %in% normFusionsDF$X.FusionName))
  if(nrow(allFusionsDFFilt) > 0)
  {
    ##Fusion genes in most of the samples
    allFusionsDFFilt$Fusion<-paste0(allFusionsDFFilt$X.FusionName,"-",allFusionsDFFilt$LeftBreakpoint,"-",allFusionsDFFilt$RightBreakpoint)
    col<-colnames(allFusionsDFFilt)
    allFusionsDFFilt$TCGAannots<-sapply(allFusionsDFFilt$annots,FUN=function(x){
      return(if(str_detect(x,"TCGA_StarF2019")) "Y" else "N" )
    })
    allFusionsDFFilt$CCLEannots<-sapply(allFusionsDFFilt$annots,FUN=function(x){
      return(if(str_detect(x,"CCLE_StarF2019")) "Y" else "N" )
    })
    uniqueFusions <- allFusionsDFFilt %>% dplyr::group_by_at(vars(Fusion)) %>% dplyr::summarize(SpanningFragCount=paste(SpanningFragCount,collapse = ","),JunctionReadCount=paste(JunctionReadCount,collapse = ","),FFPM=paste(FFPM,collapse = ","),Sample = paste(sample,collapse = ","),Number = dplyr::n(),TCGA=unique(TCGAannots),CCLE=unique(CCLEannots)) %>% arrange(desc(Number))

    #allFusionsDFFilt<-allFusionsDFFilt %>% filter(X.FusionName %in% c("RN7SL2--RF00100"))

    mergeFusionDF<-merge(uniqueFusions,allFusionsDFFilt,by="Fusion",all.x=TRUE,all.y=FALSE)
    curr_fusions<-mergeFusionDF[!duplicated(mergeFusionDF$X.FusionName),]
    left_break <- do.call(rbind,strsplit(curr_fusions$LeftBreakpoint,":"))
    right_break <- do.call(rbind,strsplit(curr_fusions$RightBreakpoint,":"))
    gene_names <- do.call(rbind,strsplit(curr_fusions$X.FusionName,"--"))
    sampleid <- (curr_fusions$Sample)
    SpanningFragCount <- curr_fusions$SpanningFragCount.x
    JunctionReadCount <- curr_fusions$JunctionReadCount.x
    ffpm <- curr_fusions$FFPM.x
    nos <- curr_fusions$Number
    TCGAAnnot<-curr_fusions$TCGAannots
    CCLEAnnot<-curr_fusions$CCLEannots

    if(nrow(curr_fusions) >1)
    {
      return_df <- as.matrix(cbind(gene_names,left_break[,1:2], right_break[,1:2],ffpm,SpanningFragCount,JunctionReadCount,nos,sampleid,TCGAAnnot,CCLEAnnot))
      colnames(return_df) <- c("Gene1","Gene2","Gene1_Chr","Gene1_Position","Gene2_Chr","Gene2_Position","FFPM","SpanningFragCount","JunctionReadCount","No. Of Samples","sample_id","TCGAAnnot","CCLEAnnot")
    }else
    {
      return_df <- t(as.matrix(c(gene_names,left_break[,1:2], right_break[,1:2],ffpm,SpanningFragCount,JunctionReadCount,nos,sampleid,TCGAAnnot,CCLEAnnot)))
      colnames(return_df) <- c("Gene1","Gene2","Gene1_Chr","Gene1_Position","Gene2_Chr","Gene2_Position","FFPM","SpanningFragCount","JunctionReadCount","No. Of Samples","sample_id","TCGAAnnot","CCLEAnnot")

    }
    links_data<-(data.frame(return_df))

    name<-str_replace_all(name,pattern = "[ /]",replacement = ".")
    write.table(links_data,file = paste0(mywd,"UniqueFusions-",name,".txt"),sep = "\t",quote = F,row.names = F)
    l[[name]]<-links_data
  }
}

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(l, paste0(mywd,"unique.fusion.genes.WONormal.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)

