library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)

# l<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/cancerSpecificGenes-TPM.rda"))

#### Inputs
mywd="/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"  ## Working directory
cnv_out_dir <- "/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"     ## Location of cnv output files
cnv_suffix <- "*driver.catalog.somatic.tsv"
sample_files <- list.files(path = cnv_out_dir,pattern = cnv_suffix,recursive = F, full.names = T)
names(sample_files) <- path_file(sample_files)#basename(dirname(basename(dirname(sample_files))))
all_cnvs_results <- lapply(sample_files, function(currfile) {
  # currfile=sample_files[[82]]
  samplename <-  unlist(str_split(path_file(currfile),pattern = ".driver"))[1]
  # print(samplename)
  curr_results <- read.table(currfile, header=T, sep="\t", comment.char = "", stringsAsFactors = F)
  if (nrow(curr_results) > 1) {
    curr_results$sample <- samplename
  } else {
    curr_results <- NA
  }
  return(curr_results)
})

all_cnvs_results <- all_cnvs_results[!is.na(all_cnvs_results)]
sNames<-unlist(lapply(names(all_cnvs_results),FUN=function(x){
  return(unlist(str_split(x,pattern = ".driver"))[1])
}))
names(all_cnvs_results)<-sNames

#Feature files
features.file <- "/Users/jaina13/myPART/WGSData/featuresWGS.txt"
features <- read.table(features.file, sep="\t", header = T, stringsAsFactors = F)
row.names(features)<-unlist(lapply(features$WGSSampleID,FUN=function(x){return(unlist(str_split(x,pattern = "ample_"))[2])}))#paste0()
features<-features[names(all_cnvs_results),]
NETSamples<-features[grepl("Neuroendocrine",features$New.Diagnosis),]
features[grepl("Neuroendocrine",features$New.Diagnosis),]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)


##CNVs analysis in individual groups
#mywd<-"/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"
# normal_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c("Control Tumor") ))#& Tissue %in% c("Adrenal")))
# normal_cnv_res <-all_cnvs_results[normal_samples]
# normCNVDF <- do.call("rbind",normal_cnv_res)
# normCNVDF$uniq<-paste0(normCNVDF$chromosome,"-",normCNVDF$chromosomeBand,"-",normCNVDF$gene,"-",normCNVDF$driver)
l<-list()
l1<-list()
names<-unique(features$New.Diagnosis)
for(name in names[names !="Control Tumor"])
{
  tumor_samples <- row.names(features %>% dplyr::filter(New.Diagnosis %in% c(name)))#,"Adrenocortical carcinoma""Gastrointestinal stromal tumor"
  tumor_cnv_res <- all_cnvs_results[tumor_samples]
  allCNVDF <- do.call("rbind",tumor_cnv_res)
  allCNVDF$uniq<-paste0(allCNVDF$chromosome,"-",allCNVDF$chromosomeBand,"-",allCNVDF$gene,"-",allCNVDF$driver)
  allCNVDFFilt<-allCNVDF #%>% dplyr::filter(!uniq %in% normCNVDF$uniq)
  uniqueCNVs <- allCNVDFFilt %>% dplyr::group_by_at(vars(uniq)) %>% dplyr::summarize(chromosome=unique(chromosome),chromosomeBand=unique(chromosomeBand),gene=unique(gene),driver=unique(driver),minCopyNumber=paste(minCopyNumber,collapse = ","),
                                                                                   maxCopyNumber=paste(maxCopyNumber,collapse = ","),
                                                                                   Sample = paste(sample,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))
  ##Extract CancerSpecificGenes
  # c<-cellSpecificGenesOutput %>% dplyr::filter(Group %in% c("Tissue-Enriched")) %>% dplyr::filter(Tissue == name)
  # c1<-merge(x=c,y=mapping,by.x="Gene",by.y="gene_id",all.x=TRUE,all.y=FALSE)
  # c2<-merge(x=uniqueCNVs,y=c1,by.x="gene",by.y="Gene",all.x=TRUE,all.y=FALSE)

  name<-str_replace_all(name,pattern = "[ /]",replacement = ".")
  #write.table(uniqueCNVs,file = paste0(mywd,name,".unique.driver.catalog.somatic.tsv"))
  l[[name]]<-uniqueCNVs %>% dplyr::select(-uniq)
  l1[[name]]<-allCNVDFFilt
  #l[[name]]<-c2 %>% dplyr::select(-uniq)
}

save(l,l1,file = paste0(mywd,"driverGenesByCancer.rda"))
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(l, paste0(mywd,"unique.driver.catalog.somatic.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)
#d<-getMAFDashboard(MAFfilePath = "/Users/jaina13/myPART/WESData/mafs/1_RT00032_W309_vs_35_RT00032_RAST_WB_0026_001.hard-filtered.maf",plotList = "burden",outputFilePath = "/Users/jaina13/myPART/WESData/mafs/",outputFileName = "1_RT00032_W309_vs_35_RT00032_RAST_WB_0026_001.html")


###Code to make karyoplot for the cancer wise samples
#https://www.biostars.org/p/257170/
features.file <- "/Users/jaina13/myPART/WGSData/featuresWGS.txt"
features <- read.table(features.file, sep="\t", header = T, stringsAsFactors = F)
features$SampleId<-unlist(lapply(features$WGSSampleID,FUN=function(x){return(unlist(str_split(x,pattern = "ample_"))[2])}))#paste0()
purpleDir<-"~/myPART/WGSData/CNVs/"
max.cnv <-8#max(unlist(lapply(purpleCNVData, function(currfile){max(currfile$copyNumber)})))
#featuresFilt<-features[features$Diagnosis == "Gastrointestinal stromal tumor",]
featuresFilt<-features[features$New.Diagnosis != "Control Tumor",]

###Purple cnv.somatic.tsv
purpleCNVDataSomatic<-apply(featuresFilt,1,function(x) {
  print(paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.cnv.somatic.tsv"))
  d<-read.table(file = paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.cnv.somatic.tsv"), sep = "\t", header = T, comment.char = "!") %>%
    mutate(chromosome = gsub("chr", "", chromosome)) %>%
    filter(!chromosome %in% c('X','Y'), bafCount > 0)
  d$chromosome<-paste0("chr",d$chromosome)
  ##There are some copuNumber that are negative
  # d$copyNumber<-abs(d$copyNumber)
  ##Making the high CNVs to 8 for plot
  # d$copyNumber<-ifelse(d$copyNumber > max.cnv, max.cnv,d$copyNumber)
  ##Calculating status based on the PURPLE manual
  purityData<-read.table(file = paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.purity.tsv"), sep = "\t", header = T, comment.char = "!")
  ploidy <- as.numeric(unlist(purityData[,"ploidy"]))
  d$status<-unlist(apply(d,1,FUN=function(x){
    minCopyNumber = as.numeric(unlist(x["copyNumber"]))
    if(minCopyNumber > (3*ploidy))
    {
      return("gain")
    }else if(minCopyNumber < 0.5)
    {
      return("loss")
    }else{
      return("neutral")
    }
  }))
  d$SampleId<-x["RTNo"]
  d$UniqVar<-paste0(d$chromosome,"-",d$start,"-",d$end)
  return(d)
})
names(purpleCNVDataSomatic)<-featuresFilt$RTNo
allCNV<-do.call(rbind,purpleCNVDataSomatic)
uniqueCNV <- allCNV %>% dplyr::group_by_at(vars(UniqVar)) %>% dplyr::summarize(ChrStartEnd=unique(UniqVar),
                                                                                     Sample = paste(`SampleId`,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))



###Purple cnv.gene.tsv
purpleCNVDataGene<-apply(featuresFilt,1,function(x) {
  print(paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.cnv.gene.tsv"))
  d<-read.table(file = paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.cnv.gene.tsv"), sep = "\t", header = T, comment.char = "!") %>%
    mutate(chromosome = gsub("chr", "", chromosome)) %>%
    filter(!chromosome %in% c('X','Y')) #minRegionMethod != "UNKNOWN" )
  d$chromosome<-paste0("chr",d$chromosome)
  ##There are some copyNumber that are negative
  #d$minCopyNumber<-abs(d$minCopyNumber)

  ##Making the high CNVs to 8 for plot
  #d$minCopyNumber<-ifelse(d$minCopyNumber > max.cnv, max.cnv,d$minCopyNumber)
  ##Calculating status based on the PURPLE manual
  purityData<-read.table(file = paste0(purpleDir, "/", x["WGSSampleID"],"/purple/",x["SampleId"], ".purple.purity.tsv"), sep = "\t", header = T, comment.char = "!")
  ploidy <- as.numeric(unlist(purityData[,"ploidy"]))
  d$status<-unlist(apply(d,1,FUN=function(x){
    minCopyNumber = as.numeric(unlist(x["minCopyNumber"]))
    if(minCopyNumber > (3*ploidy))
    {
      return("gain")
    }else if(minCopyNumber < 0.5)
    {
      return("loss")
    }else{
      return("neutral")
    }
  }))
  d$SampleId<-x["RTNo"]
  #d$UniqVar<-paste0(d$chromosome,"-",d$start,"-",d$end)
  return(d)
})

names(purpleCNVDataGene)<-featuresFilt$RTNo
allCNV<-do.call(rbind,purpleCNVDataGene)
uniqueCNV <- allCNV %>% dplyr::group_by_at(vars(UniqVar)) %>% dplyr::summarize(ChrStartEnd=unique(UniqVar),
                                                                               Sample = paste(`SampleId`,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))
##Make the matrix for CNVs across all samples
copyNumberList<-lapply(purpleCNVDataGene,FUN = function(x){
  row.names(x)<-x$gene
  return(x[,5,drop=FALSE])
})
copyNumber<-do.call(cbind,copyNumberList)
colnames(copyNumber)<-names(purpleCNVDataGene)
copyNumberFilt<-copyNumber[unique(allCNV[allCNV$status != "neutral",]$gene),]
write.table(copyNumberFilt,paste0(mywd,"../minCopyNumberMatrixDriverGenes.txt"),quote = FALSE,sep = "\t")


#https://github.com/hartwigmedical/hmftools/issues/88
###Code to group the genes based on Purple Results
autosomes <- paste("chr", 1:22, sep = '')
samplesCancerPurity<<-c()
copyTable <- do.call(base::rbind, apply(featuresFilt,1, function(sampleID)
{
  #print(sampleID)
  #print(is_file(paste0(purpleDir, "/", sampleID["WGSSampleID"],"/purple/",sampleID["SampleId"], ".purple.cnv.gene.tsv")))
  geneCopyNumbers <- read.delim((paste0(purpleDir, "/", sampleID["WGSSampleID"],"/purple/",sampleID["SampleId"], ".purple.cnv.gene.tsv")))
  #print(paste0(purpleDir, "/", sampleID["WGSSampleID"],"/purple/",sampleID["SampleId"], ".purple.purity.tsv"))
  samplePurityPloidy <- read.delim((paste0(purpleDir, "/", sampleID["WGSSampleID"],"/purple/",sampleID["SampleId"], ".purple.purity.tsv")))
  sampleOverallPloidy <- round(samplePurityPloidy[["ploidy"]])
  samplesCancerPurity <<- c(samplesCancerPurity,samplePurityPloidy[["purity"]])
  #print(class(samplePurityPloidy[["purity"]]))
  sampleCopyTable <- data.frame(geneCopyNumbers[, "chromosome"], geneCopyNumbers[, "gene"], sampleID["SampleId"],#sampleID,
                                geneCopyNumbers[, "maxCopyNumber"], stringsAsFactors = FALSE)
  sampleCopyTable[, 4] <- round(sampleCopyTable[, 4])

  copyGroups <- NA
  if(samplePurityPloidy[, "gender"] == "FEMALE")
  {
    if(sampleOverallPloidy <= 2) # Diploid
    {
      copyGroups[sampleCopyTable[, 4] == 0 & (sampleCopyTable[, 1] %in% c(autosomes, "chrX"))] <- "Deletion" # No copies left.
      copyGroups[sampleCopyTable[, 4] >= 5 & sampleCopyTable[, 1] %in% c(autosomes, "chrX")] <- "Amplification"
    } else { # Genome doubling
      copyGroups[sampleOverallPloidy - sampleCopyTable[, 4] >= 2 & (sampleCopyTable[, 1] %in% c(autosomes, "chrX"))] <- "Deletion"
      copyGroups[sampleCopyTable[, 4] >= 8 & (sampleCopyTable[, 1] %in% c(autosomes, "chrX"))] <- "Amplification"
    }
  } else { # Male
    if(sampleOverallPloidy <= 2) # Diploid
    {
      copyGroups[sampleCopyTable[, 4] == 0] <- "Deletion" # No copies left.
      copyGroups[sampleCopyTable[, 4] >= 2 & (sampleCopyTable[, 1] %in% c("chrX", "chrY"))] <- "Amplification"
      copyGroups[sampleCopyTable[, 4] >= 5 & (sampleCopyTable[, 1] %in% autosomes)] <- "Amplification"
    } else { # Genome doubling
      copyGroups[sampleOverallPloidy - sampleCopyTable[, 4] >= 3 & (sampleCopyTable[, 1] %in% autosomes)] <- "Deletion"
      copyGroups[(sampleOverallPloidy / 2) - sampleCopyTable[, 4] >= 1 & (sampleCopyTable[, 1] %in% c("chrX", "chrY"))] <- "Deletion"
      copyGroups[sampleCopyTable[, 4] >= 8 & (sampleCopyTable[, 1] %in% autosomes)] <- "Amplification"
      copyGroups[sampleCopyTable[, 4] >= 3 & (sampleCopyTable[, 1] %in% c("chrX", "chrY"))] <- "Amplification"
    }
  }

  sampleCopyTable[, 4] <- copyGroups
  sampleCopyTable <- sampleCopyTable[!is.na(sampleCopyTable[, 4]), ]
  sampleCopyTable[, 2:4]
}))

cnvMatrix<-data.frame(copyTable %>% tidyr::pivot_wider(names_from = sampleID..SampleId..,values_from=geneCopyNumbers....maxCopyNumber..))
row.names(cnvMatrix)<-cnvMatrix[,1]
cnvMatrix<-cnvMatrix[,2:ncol(cnvMatrix)]
colnames(cnvMatrix) <- gsub("^X","",colnames(cnvMatrix))
CNV_Threshold<-0.5
#keep.rows <- apply(cnvMatrix,MARGIN = 1, function(x){if(sum(is.na(x)) >= length(x)*(1-CNV_Threshold)){return(FALSE)}else{return(TRUE)}})
keep.rows <- apply(cnvMatrix,MARGIN = 1, function(x){if(sum(is.na(x)) >= (length(x)-5)){return(FALSE)}else{return(TRUE)}})
cnvMatrixFilt <- cnvMatrix[keep.rows,]
hubGenes<-analyzePPINetwork(row.names(cnvMatrixFilt),paste0(purpleDir,"cnvGenes.PPI.900.png"),ppiEdgeThreshold = 900)
write.table(hubGenes,paste0(purpleDir,"/cnvGenes.PPI.900.hubGenes.txt"),quote = F,row.names = F,sep = "\t")

# bins <- toGRanges(do.call(rbind,purpleCNVData)[,c(1,2,3)])
# nsamples <- length(purpleCNVData)

# geneRTNo <- allCNV %>% dplyr::group_by_at(vars(gene)) %>% dplyr::summarize(copyNumber=paste(`minCopyNumber`,collapse = ","),Sample = paste(`SampleId`,collapse = ","))
# allGenes<-unique(allCNV$gene)
# cnvCountMatrix <- data.frame(matrix(0,ncol = nrow(featuresFilt),nrow=length(allGenes)))
# row.names(cnvCountMatrix) <- allGenes
# colnames(cnvCountMatrix)<-featuresFilt$RTNo
# cnvMatrixList <- apply(geneRTNo,1,FUN = function(x){
#   samples<-unique(unlist(stringr::str_split(x[3],pattern = ",")))
#   copyNumber<-unique(unlist(stringr::str_split(x[2],pattern = ",")))
#   if(length(copyNumber) != length(samples))
#   {
#     print(x[1])
#   }
#   names(copyNumber)<-samples
#   cnvCountMatrix[x[1],samples]<-copyNumber[samples]
#   return(cnvCountMatrix[x[1],])
# })
# binaryCountDF<-do.call(rbind,cnvMatrixList)

##Plotting CNVs in terms of oncoplot
mywd="/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"
o<-load(file = paste0(mywd,"driverGenesByCancer.rda"))
CNVs<-rbind(l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,l1$`Neuroendocrine-lung`)
#CNVs<-l1$Adrenocortical.carcinoma
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
#featuresCNV<-features[features$RTNo %in% unique(custom.cn.data$Sample_name),]
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
wesSNVs<-l1$Neuroendocrine.tumor
wesSNVs@data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#wesSNVs@clinical.data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@clinical.data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#head(custom.cn.data)
library(openxlsx)
features<-read.xlsx("~/myPART/NGS clinical data for Ashish 4-26-21v2-no_name.xlsx")
featuresFilt <- features[features$Subject.Code %in% c(unique(custom.cn.data$Sample_name),unique(wesSNVs@data$Tumor_Sample_Barcode)),]
featuresFilt$Tumor_Sample_Barcode<-featuresFilt$Subject.Code
featuresFilt$PrimaryLocation<-featuresFilt$Location.of.Primary
featuresFilt$Tissue<-featuresFilt$`anatomical.location.-.MyPART.Tumor.Pathologies`
featuresFilt$SampleType<-featuresFilt$`sample.type.-.MyPART.Tumor.Pathologies`
featuresFilt$Age<-featuresFilt$`age.at.sample.collection.-.Path.Rpt`

maf.plus.cn = read.maf(maf = wesSNVs@data,cnTable = custom.cn.data,clinicalData = featuresFilt,verbose = FALSE)
summary<-maf.plus.cn@gene.summary
summary$Altertotal<-summary$MutatedSamples+summary$AlteredSamples
summary$Altertotal[summary$MutatedSamples != 0 & summary$CNV_total == 0]<-summary$MutatedSamples[summary$MutatedSamples != 0 & summary$CNV_total == 0]
summary<-summary[!(summary$Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")) & summary$Altertotal > 1,]
genes<-summary[order(Altertotal,decreasing = T),]$Hugo_Symbol[1:50]

pdf(paste0(mywd,"SNVs-CNVs-NETs.pdf"),width = 10,height = 12)
oncoplot(maf = maf.plus.cn,genes = genes,legendFontSize = 1.2,clinicalFeatures = c("PrimaryLocation","Tissue","SampleType"),gene_mar = 10)
dev.off()

uniqueCNVs <- custom.cn.data %>% dplyr::group_by_at(vars(Gene)) %>%
  dplyr::summarize(Gene=unique(Gene),CN=unique(CN),Sample = paste(Sample_name,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))
