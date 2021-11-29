###########Figure 3A##################
mywd="/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"
o<-load(file = paste0(mywd,"driverGenesByCancer.rda"))
#CNVs<-rbind(l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,l1$`Neuroendocrine-lung`)
CNVs<-rbind(l1$Adrenocortical.carcinoma,l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,
            l1$`Neuroendocrine-lung`,l1$Pheochromocytoma,l1$Medullary.thyroid.carcinoma,l1$Papillary.thyroid.carcinoma,l1$Anaplastic.thyroid.carcinoma)
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
#featuresCNV<-features[features$RTNo %in% unique(custom.cn.data$Sample_name),]

###CNVs from Sequenza from WES data
o<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/CNV/sequenza_out/sequenzaRes-Unfilt.rda"))
CNVs<-rbind(cn.list.amp$`11_RT00050_W305`,cn.list.amp$`29_RT00142_X646`,cn.list.amp$`19_RT00108_X005`)
CNVs$RTNo<-unlist(lapply(as.character(CNVs$Sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$Type,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "Deletion"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "Amplification"]<-"Amp"

# o2<-load(file = "/Users/jaina13/myPART/WGSData/tumor-only-somatic-mafs/filterMAFResults.rda")
# wgsSNVs<-l1$Adrenocortical.carcinoma
# wgsSNVs@data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wgsSNVs@data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[2]}))
# d<-merge(maf.plus.cn@variants.per.sample,sample_annotation_data,by="Tumor_Sample_Barcode")

o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/merged_somatic_variants/maf-new/SNVsResults-consensus.rda"))
#o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
#wesSNVs<-l1$Adrenocortical.carcinoma
wesSNVs<-maftools::merge_mafs(list(l1$Adrenocortical.carcinoma,l1$Neuroendocrine.tumor,l1$Pheochromocytoma,l1$Medullary.thyroid.carcinoma,l1$Anaplastic.thyroid.cancer))
wesSNVs@data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[2]}))
#wesSNVs<-read.maf(wesSNVs@data)
#wesSNVs@clinical.data$Tumor_Sample_Barcode<-unlist(lapply(as.character(wesSNVs@clinical.data$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[3]}))
#head(custom.cn.data)
library(openxlsx)
features<-openxlsx::read.xlsx("~/myPART/NGS clinical data for Ashish 4-26-21v2-no_name.xlsx")
featuresFilt <- features[features$Subject.Code %in% c(unique(custom.cn.data$Sample_name),unique(wesSNVs@data$Tumor_Sample_Barcode)),]
#featuresFilt <- features[features$Subject.Code %in% unique(custom.cn.data$Sample_name),]
featuresFilt$Tumor_Sample_Barcode<-featuresFilt$Subject.Code
featuresFilt$PrimaryLocation<-tools::toTitleCase(featuresFilt$Location.of.Primary)
featuresFilt$Tissue<-tools::toTitleCase(featuresFilt$`anatomical.location.-.MyPART.Tumor.Pathologies`)
featuresFilt$SampleType<-tools::toTitleCase(featuresFilt$`sample.type.-.MyPART.Tumor.Pathologies`)
featuresFilt$Age<-featuresFilt$`age.at.sample.collection.-.Path.Rpt`

# featuresFilt %>% dplyr::group_by_at(vars(`cancer/tumor.name.-.MyPART.Tumor.Pathologies`)) %>% dplyr::summarize(gene=unique(`Hugo Symbol`),chromosome=unique(Chromosome),start=unique(`Start Position`),end=unique(`End Position`),
#                                                                                  VariantType =paste(unique(`Variant Type`),collapse = ","),
#                                                                                  TranscriptChange=paste(unique(`Transcript Change`),collapse = ","),ProteinChange=paste(unique(`Protein Change`),collapse = ","),
#                                                                                  ExistingAnnotation=paste(unique(`Existing Annotation`),collapse = ","),EffectPredictionSIFT=paste(unique(`Effect Prediction - SIFT`),collapse = ","),
#                                                                                  EffectPredictionPolyPhen=paste(unique(`Effect Prediction - PolyPhen`),collapse = ","),ClinVar=paste(unique(`Known Effects ClinVar`),collapse = ","),
#                                                                                  Sample = paste(`Sample ID`,collapse = ","),NumberOfSamples = dplyr::n()) %>% arrange(desc(NumberOfSamples))

lapply(unique(featuresFilt$New.Diagnosis),function(x){
  f <- featuresFilt[featuresFilt$New.Diagnosis == x,]
  print(x)
  print(table(f$Race))
  #print(sd(as.numeric(f$Age),na.rm=TRUE))
})

#maf.plus.cn = read.maf(maf = wesSNVs@data[ wesSNVs@data$Hugo_Symbol== "HTR1D",],cnTable = custom.cn.data,clinicalData = featuresFilt,verbose = FALSE)
#maf.plus.cn = read.maf(maf = wesSNVs@data,clinicalData = featuresFilt,verbose = FALSE)
maf.plus.cn = read.maf(maf = wesSNVs@data,cnTable = custom.cn.data,clinicalData = featuresFilt,verbose = FALSE)

# summary<-maf.plus.cn@gene.summary
# summary$Altertotal<-summary$MutatedSamples+summary$CNV_total
# summary$Altertotal[summary$MutatedSamples != 0 & summary$CNV_total == 0]<-summary$MutatedSamples[summary$MutatedSamples != 0 & summary$CNV_total == 0]
# #summary<-summary[!(summary$Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")) & summary$Altertotal > 1,]
# summary<-summary[!(summary$Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")),]
#genes<-summary[order(Altertotal,decreasing = T),]$Hugo_Symbol[1:50]
#variant.type.summary <- maf.plus.cn@variant.type.summary
#variant.type.summary$total <- variant.type.summary$total + variant.type.summary$CNV
#maf.plus.cn@gene.summary <- summary
#maf.plus.cn@variant.type.summary <- variant.type.summary

# source("~/maftools/R/oncomatrix.R")
# source("~/maftools/R/oncoplot.R")
# pdf(paste0(mywd,"SNVs-CNVs-Pheocromocytoma.pdf"),width = 10,height = 12)
# maftools::oncoplot(maf = maf.plus.cn,genes = genes,legend_height = 10,legendFontSize = 2,annotationFontSize = 2,clinicalFeatures = c("PrimaryLocation","Tissue","SampleType"),gene_mar = 10,barcode_mar = 10)
# dev.off()

source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/Figure3A-helper.R")
sample_annotation_data<-maf.plus.cn@clinical.data[,c("Tumor_Sample_Barcode","cancer/tumor.name.-.MyPART.Tumor.Pathologies","SampleType","Tissue","Age","TMB.(Mut/Mb).-.MyPART.TSO.500s")]
colnames(sample_annotation_data)<-c("Tumor_Sample_Barcode","Diagnosis","SampleType","Tissue","Age","TMB_TSO500")
sample_annotation_data$Age <-as.numeric(sample_annotation_data$Age)
sample_annotation_data$TMB_TSO500<-as.numeric(sample_annotation_data$TMB_TSO500)
sample_annotation_colors <- list(Diagnosis=diagnosisColorDNASeq[sample_annotation_data$Diagnosis],SampleType=tumorSiteColor[sample_annotation_data$SampleType],Tissue=tissueColorDNASeq[sample_annotation_data$Tissue],Age=getAgeColor(sample_annotation_data$Age),TMB_TSO500 = getTMB500Color((sample_annotation_data$TMB_TSO500)))#get_clinical_colors(sample_annotation_data)
# sample_annotation_data<-maf.plus.cn@clinical.data[,c("Tumor_Sample_Barcode","cancer/tumor.name.-.MyPART.Tumor.Pathologies","SampleType","Tissue","Age")]
# colnames(sample_annotation_data)<-c("Tumor_Sample_Barcode","Diagnosis","SampleType","Tissue","Age")
# sample_annotation_data$Age <-as.numeric(sample_annotation_data$Age)
# sample_annotation_colors <- list(Diagnosis=diagnosisColorDNASeq[sample_annotation_data$Diagnosis],SampleType=tumorSiteColor[sample_annotation_data$SampleType],Tissue=tissueColorDNASeq[sample_annotation_data$Tissue],Age=getAgeColor(sample_annotation_data$Age))#get_clinical_colors(sample_annotation_data)

g<-make_oncoplot(maf.plus.cn,show_sample_names = FALSE,clin_data = sample_annotation_data,clin_data_colors = sample_annotation_colors,ngene_max = 50)
#pdf(paste0(mywd,"CNVs-ALL-Endocrines.pdf"),width = 14,height = 12)
pdf(paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/merged_somatic_variants/maf-new/","SNVs-consensus-SequenzaCNV.pdf"),width = 14,height = 12)
#ComplexHeatmap::draw(g, show_annotation_legend = TRUE)
draw(onco_base_default, show_annotation_legend = TRUE)
dev.off()



make_oncoplot <- function(maf.filtered, cohort_freq_thresh = 0.01,ngene_max=25, auto_adjust_threshold=T,
                          oncomat_only=F, show_sample_names=NULL,
                          clin_data=NULL, clin_data_colors=NULL,
                          savename=NULL) {

  maf.filtered = maf.plus.cn
  clin_data=sample_annotation_data
  clin_data_colors=sample_annotation_colors
  cohort_freq_thresh = 0.01
  oncomat_only=F
  show_sample_names=FALSE
  auto_adjust_threshold=T
  ngene_max=50

  require(ComplexHeatmap)

  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$AlteredSamples/as.numeric(maf.filtered@summary$summary[3])),
                         mutation_count=maf.filtered@gene.summary$total + maf.filtered@gene.summary$CNV_total,
                         #mutation_count=maf.filtered@gene.summary$total,#For SNV only
                         #mutation_count=maf.filtered@gene.summary$CNV_total, #For CNV only
                         stringsAsFactors = F)

  #frac_mut <- frac_mut %>% dplyr::filter(!(Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C","HTR1D")))
  frac_mut <- frac_mut %>% dplyr::filter(!(Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")) & (mutation_count > 1))

  target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
  if (auto_adjust_threshold) {
    cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
  }
  ### Select genes based on the frequency threshold
  #frac_mut <- frac_mut[order(frac_mut$frac_mut,frac_mut$mutation_count,decreasing = T),]
  frac_mut <- frac_mut[order(frac_mut$frac_mut,decreasing = T),]
  freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
  freq_genes <- freq_genes[1:min(ngene_max,length(freq_genes))]
  if (length(freq_genes) == 0) {
    stop("No genes to plot; change the frequency threshold to include more genes.")
  }
  if (length(freq_genes) > 100) {
    target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))],2)
    # stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    warning(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    # return(NA)
  }
  gene_list <- list(freq_genes)
  reasons <- paste0("Cohort Freq > ",round(cohort_freq_thresh,digits = 3))

  ### Collect genes to plot
  genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
  for (i in 1:length(gene_list)) {
    if (is.na(gene_list[[i]][1])) {
      next
    }
    genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                data.frame(Hugo_Symbol=gene_list[[i]],
                                           reason=reasons[i]))
  }
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])

  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]
  # browser()
  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=genes_for_oncoplot$reason
  split_colors <- rainbow(length(unique(split_idx)))
  # names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
  names(split_colors) <- unique(split_idx)
  split_colors <- list(Reason=split_colors)

  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  #print(genes_for_oncoplot$Hugo_Symbol)
  oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
  onco_genes <- rownames(oncomat)

  if (oncomat_only) {
    return(oncomat)
  }
  oncomat.plot <- oncomat

  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
  }

  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)
  # browser()
  if (length(oncomat.plot) < 1) {
    stop("No samples to plot.")
  }
  ### Column labels get cluttered if too many samples
  if (!is.logical(show_sample_names)) {
    show_sample_names=T
    if (ncol(oncomat.plot) > 20) {
      show_sample_names=F
    }
  }

  myanno=NULL
  if (!is.null(clin_data)) {
    # browser()
    anno_data <- data.frame(clin_data[match(colnames(oncomat.plot), clin_data$Tumor_Sample_Barcode),],stringsAsFactors = F)
    row.names(anno_data) <- anno_data$Tumor_Sample_Barcode
    anno_data <- anno_data[,!colnames(anno_data) %in% "Tumor_Sample_Barcode", drop=F]
    if (ncol(anno_data) > 0) {
      ###Make changes in the text for the annotation
      myanno <- HeatmapAnnotation(df=anno_data,col = clin_data_colors,
                                  annotation_name_gp =  gpar(fontsize = 14,fontface = 2),
                                  annotation_legend_param=list(labels_gp = gpar(fontsize = 14),title_gp = gpar(fontsize = 16, fontface = 2),
                                                               nrow = 5,
                                                               legend_direction = "vertical"))
    }
  }

  ## Show total burden for top annotation
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  #variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[,!names(variant_type_data)%in%c("Tumor Sample Barcode","total","CNV total")]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  ###Make changes in the text/data for the top mutation histogram
  top_ha = HeatmapAnnotation("Total\nMutations" = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F),
                             annotation_name_side = "left",annotation_name_rot=90,annotation_name_gp = gpar(fontsize = 12),height = unit(4, "cm"))

  # browser()

  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(fontsize = 12)), show_annotation_name=F)
  # print(oncomat.plot)
  ### Make the oncoplot
  # alter_fun$Deletion<-function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
  #             gp = gpar(fill = mutation_colors["Del"], col = NA))
  # }
  # alter_fun$Amplification<-function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
  #             gp = gpar(fill = mutation_colors["Amp"], col = NA))
  # }
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 column_order = row.names(anno_data[order(anno_data$Diagnosis),]),
                                 show_pct = F,
                                 row_split=split_idx,
                                 row_title = NULL,
                                 bottom_annotation = myanno,
                                 top_annotation = top_ha,
                                 left_annotation = left_ha,
                                 show_column_names = show_sample_names,
                                 alter_fun_is_vectorized = T,
                                 row_names_gp = gpar(fontsize = 14),
                                 heatmap_legend_param = list(title = "Alterations",title_gp = gpar(fontsize = 16, fontface = 2),labels_gp = gpar(fontsize = 14)))#,

  draw(onco_base_default, show_annotation_legend = TRUE)

  if ( ! is.null(savename) ) {
    # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
    # browser()
    anno_height=ifelse(!is.null(clin_data), min(c(4, 0.5*ncol(clin_data))), 0)
    onco_height=max(round(0.1*nrow(oncomat.plot),0),4) + anno_height
    # onco_width=onco_height*0.75 + anno_height*1.2
    onco_width=ncol(oncomat.plot)*0.01 + anno_height*1.2
    pdf(file = savename,height=onco_height,width=onco_width)
    draw(onco_base_default)
    dev.off()
  }

  ### Return the oncoplot (if function is pointed to a variable)
  invisible(onco_base_default)
}

# o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
# allMAFDF<-do.call(rbind,lapply(names(l), function(n){o<-l[[n]];o$Diagnosis<-n;return(o)}))
# allMAFDF$UniqVar<-paste0(allMAFDF$chromosome,"-",allMAFDF$start,"-",allMAFDF$end,"-",allMAFDF$TranscriptChange,"-",allMAFDF$VariantType,"-",allMAFDF$Diagnosis)
# mafPipelinerAllFilt<-allMAFDF
getMAFDF<-function(l)
{
  allMAFDF<-do.call(rbind,lapply(names(l), function(n){o<-l[[n]];o$Diagnosis<-n;return(o)}))
  allMAFDF$UniqVar<-paste0(allMAFDF$chromosome,"-",allMAFDF$start,"-",allMAFDF$end,"-",allMAFDF$TranscriptChange,"-",allMAFDF$VariantType,"-",allMAFDF$Diagnosis)
  return(allMAFDF)
}

getMAFDFSample<-function(l1,RTCol=2)
{
  allMAFDF<-do.call(rbind,lapply(names(l1), function(n){o<-data.frame(l1[[n]]@data);o$Tumor_Sample_Barcode<-unlist(lapply(as.character(o$Tumor_Sample_Barcode), function(x){unlist(str_split(x,pattern = "_"))[RTCol]}));return(o)}))
  allMAFDF$UniqVar<-paste0(allMAFDF$Chromosome,"-",allMAFDF$Start_Position,"-",allMAFDF$End_Position,"-",allMAFDF$Transcript_Change,"-",allMAFDF$Variant_Type,"-",allMAFDF$Tumor_Sample_Barcode)
  return(allMAFDF)
}

o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new-FilterNCall2.rda"))
mafPipelinerN2Filt<-getMAFDFSample(l1,RTCol = 3)
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new-FilterMutect2.rda"))
mafPipelinerMutect2Filt<-getMAFDFSample(l1,RTCol = 3)


o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/sobdetector/maf/SNVsResults-consensus.rda"))
mafFFPEN2Filt<-getMAFDFSample(l1)
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/sobdetector/maf/SNVsResults-mutect2.rda"))
mafFFPEMutect2Filt<-getMAFDFSample(l1)


o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/merged_somatic_variants/maf/SNVsResults-consensous.rda"))
mafN2Filt<-getMAFDFSample(l1)
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/new-exome-pipeline-results/merged_somatic_variants/maf/SNVsResults-mutect2.rda"))
mafMutect2Filt<-getMAFDFSample(l1)


# install.packages("ggVennDiagram")
library(ggVennDiagram)
x <- list(mafPipelinerN2Filt= mafPipelinerN2Filt$UniqVar, mafPipelinerMutect2Filt = mafPipelinerMutect2Filt$UniqVar, mafN2Filt = mafN2Filt$UniqVar, mafMutect2Filt = mafMutect2Filt$UniqVar)
x <- list(mafPipelinerMutect2Filt = mafPipelinerMutect2Filt$UniqVar, mafMutect2Filt = mafMutect2Filt$UniqVar)
x <- list(mafFFPEN2Filt = mafFFPEN2Filt$UniqVar, mafN2Filt = mafN2Filt$UniqVar)
ggVennDiagram(x)

