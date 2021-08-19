library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)
library(maftools)

#### Inputs
mywd="/Users/jaina13/myPART/TSO500Data-POB-ACC/"  ## Working directory
#maf_out_dir <- "/Users/jaina13/myPART/WGSData/tumor-only-somatic-mafs/"     ## Location of cnv output files
snv_suffix <- "*variants.txt"
sample_files <- list.files(path = mywd,pattern = snv_suffix,recursive = F, full.names = T)
names(sample_files) <- path_file(sample_files)#basename(dirname(basename(dirname(sample_files))))
#fileNames<-unlist(lapply(sample_files, function(currfile){unlist(str_split(path_file(currfile),pattern = ".hard-filtered"))[1]}))

snvs_results_list <- lapply(sample_files, function(currfile) {
  #samplename <-  unlist(str_split(path_file(currfile),pattern = ".hard-filtered"))[1]
  #curr_results <- read.maf(paste0(maf_out_dir,currfile))
  curr_result<-read.table(currfile,sep = "\t",header = T)
  # d<-generateVariantTableFromMAFData(curr_results)
  return(curr_result[,!colnames(curr_result) %in% c("QCI.Assessment","QCI.Actionability","QCI.Nof1.Actionability")])
  # return(d)
})
snvs_results<-do.call(rbind,snvs_results_list)

##Oncogenomics ID Mapping
labMatrixMapping<-read.table(paste0(mywd,"MyPARTTSO500-OncogenomicsMapping-LabMatrix.csv"),sep = ",",header = T)
oncoMetaData<-read.table(paste0(mywd,"26001_POB.meta.txt"),sep = "\t",header = T)
mapping<-apply(oncoMetaData,1,FUN=function(x)
  {
    biomaterialIds<-unlist(stringr::str_split(x["Source.Biomaterial.ID"],pattern = ","))[1]
    #print(biomaterialIds)
    m<-labMatrixMapping[grepl(biomaterialIds,labMatrixMapping$DNAID...MyPART.TSO.500s),]
    #print(class(m))
    if(nrow(m)>0)
    {
      m$oncoPatientId<-x["Patient.ID"]
      m$Diagnosis<-x["Diagnosis"]
      return(m)
    }
    return(NULL)
})
mappingDF<-do.call(rbind,mapping[lengths(mapping) != 0])

# f<-read.xlsx("~/myPART/NGS clinical data for Ashish 4-26-21v2-no_name.xlsx")
# RTCancerMapping<-f[,c("Subject.Code","cancer/tumor.name.-.MyPART.Tumor.Pathologies")]
# m<-merge(mappingDF,RTCancerMapping,by.x="Subject.Code....Subject.Data",by.y="Subject.Code",all.x=TRUE,all.y=FALSE)

#maf<-read.maf(snvs_results)
snvsTable<-generateVariantTableFromMAFData(snvs_results)
snvsTable<-snvsTable %>% dplyr::filter(!(`Variant Classification` %in% c("Benign","Likely benign","")))
snvsTable$uniqVar<-paste0(snvsTable$Chromosome,"-",snvsTable$`Start Position`,"-",snvsTable$`End Position`,"-",snvsTable$`Transcript Change`,"-",snvsTable$`Variant Type`)
uniqueMAF <- snvsTable %>% dplyr::group_by_at(vars(uniqVar)) %>% dplyr::summarize(gene=unique(`Hugo Symbol`),chromosome=unique(Chromosome),start=unique(`Start Position`),end=unique(`End Position`),
                                                                                 VariantType =paste(unique(`Variant Type`),collapse = ","),
                                                                                 TranscriptChange=paste(unique(`Transcript Change`),collapse = ","),ProteinChange=paste(unique(`Protein Change`),collapse = ","),
                                                                                 VariantClassification=paste(unique(`Variant Classification`),collapse = ","),ClinVar=paste(unique(`Known Effects ClinVar`),collapse = ","),
                                                                                 Sample = paste(`Sample ID`,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(uniqueMAF, paste0(mywd,"unique.tso500.somatic.mutations.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)


###CNVs from the TSO500 data
cnv_suffix <- "*CNVs.csv"
sample_files <- list.files(path = mywd,pattern = cnv_suffix,recursive = F, full.names = T)
names(sample_files) <- path_file(sample_files)#basename(dirname(basename(dirname(sample_files))))
cnvs_results_list <- lapply(sample_files, function(currfile) {
  curr_result<-read.table(currfile,sep = ",",header = T)
  return(curr_result)
})
cnvs_results<-do.call(rbind,cnvs_results_list)


generateVariantTableFromMAFData <- function(maf, use_syn=FALSE, extra_cols=c()) {

  output_data<-maf
  if (all(c("Tumor_Seq_Allele1","Tumor_Seq_Allele2") %in% colnames(output_data))) {
    output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
  }
  if (all(c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2") %in% colnames(output_data))) {
    output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")
  }

  # if (! "tumor_freq" %in% colnames(output_data)) {
  #   # browser()
  #   if (all(c("t_depth","t_alt_count")%in% colnames(output_data))) {
  #     output_data$tumor_freq <- as.numeric(as.character(output_data$t_alt_count))/as.numeric(as.character(output_data$t_depth))
  #   }
  # }
  cols_for_table <- c("Hugo Symbol" = "Gene",
                      "Sample ID" = "Sample.ID",
                      "Variant Classification"="Intervar",
                      "Chromosome"="Chr","Start Position" ="Start","End Position"="End",
                      #"Reference Allele"="Ref",
                      "Variant Type"="Exonic.function",
                      "Tumor Genotype"="Alt",
                      "Normal Genotype"="Ref",
                      "Known Effects ClinVar"="Clinvar",
                      "Transcript Change"="Transcript",
                      "Protein Change"="AAChange",
                      "Tumor Alt Frequency"="VAF"
  )
  # browser()
  cols_for_table <- c(cols_for_table, extra_cols)
  # variant_info <- as.data.frame(output_data)[,cols_for_table]
  # norm_info_cols <- grep("^n_",cols_for_table, value=T)
  # # mydat <- apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)})
  # if (sum(rowSums(apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)}), na.rm=T), na.rm = T)==0) {
  #   cols_for_table <- cols_for_table[!cols_for_table %in% norm_info_cols]
  # }
  output_cols <- colnames(output_data)[match(cols_for_table, colnames(output_data), nomatch=0)]
  not_output <- cols_for_table[!cols_for_table %in% output_cols]
  if (length(not_output) > 0) {
    warning(paste0("Not outputting these columns: ", paste(not_output, collapse=", ")))
  }
  variant_info <- as.data.frame(output_data)[,output_cols]
  colnames(variant_info) <- names(cols_for_table)[match(colnames(variant_info),cols_for_table)]
  return(variant_info)
}
###PTPRN2 genes in NET-SI samples locations
samples<-"18_RT00086_X212,42_RT00085_X007,52_RT00110_W998,59_RT00185_Y172,64_RT00182_Y095,67_RT00178_X967,8_RT00041_W042_resent"
purpleDir<-"~/myPART/WGSData/CNVs/"
cnvsLocations<-lapply(unlist(str_split(samples,pattern = ",")),FUN=function(x){
  #driverGeneFile<-read.table(paste0(purpleDir, "/Sample_", x,"/purple/",x, ".driver.catalog.somatic.tsv"))
  print(paste0(purpleDir, "/Sample_", x,"/purple/",x, ".purple.cnv.gene.tsv"))
  geneCopyNumbers <- read.table((paste0(purpleDir, "/Sample_", x,"/purple/",x, ".purple.cnv.gene.tsv")),sep = "\t",header = T)
  minStartRegion<-geneCopyNumbers[geneCopyNumbers$gene == "PTPRN2",]$minRegionStart
  minStopRegion<-geneCopyNumbers[geneCopyNumbers$gene == "PTPRN2",]$minRegionEnd
  print(paste0(purpleDir, "/Sample_", x,"/purple/",x, ".purple.cnv.somatic"))
  somaticCopyNumbers <- read.table((paste0(purpleDir, "/Sample_", x,"/purple/",x, ".purple.cnv.somatic.tsv")),sep = "\t",header = T)
  somaticCopyNumbers$SampleNumber<-x
  return(somaticCopyNumbers[somaticCopyNumbers$start == minStartRegion,])
})
cnvsLocationsDF<-do.call(rbind,cnvsLocations)


