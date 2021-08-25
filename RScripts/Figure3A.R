###########Figure 3A##################
mywd="/Users/jaina13/myPART/WGSData/CNVs/DriverGenes/"
o<-load(file = paste0(mywd,"driverGenesByCancer.rda"))
#CNVs<-rbind(l1$`Neuroendocrine-small.intestine`,l1$`Neuroendocrine-pancreas`,l1$`Neuroendocrine-lung`)
CNVs<-l1$Pheochromocytoma
CNVs$RTNo<-unlist(lapply(as.character(CNVs$sample), function(x){unlist(str_split(x,pattern = "_"))[2]}))
custom.cn.data = data.frame(Gene = CNVs$gene,Sample_name = CNVs$RTNo,CN = CNVs$driver,stringsAsFactors = FALSE)
custom.cn.data$CN[custom.cn.data$CN == "DEL"]<-"Del"
custom.cn.data$CN[custom.cn.data$CN == "AMP"]<-"Amp"
custom.cn.data$CN[custom.cn.data$CN == "PARTIAL_AMP"]<-"Amp"
#featuresCNV<-features[features$RTNo %in% unique(custom.cn.data$Sample_name),]
o1<-load(file = paste0("/Users/jaina13/myPART/WESData/Pipeliner_somaticpairs/merged_somatic/SNVsResults-new.rda"))
wesSNVs<-l1$Pheochromocytoma
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
#summary<-summary[!(summary$Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")) & summary$Altertotal > 1,]
summary<-summary[!(summary$Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")),]
genes<-summary[order(Altertotal,decreasing = T),]$Hugo_Symbol[1:50]


# source("~/maftools/R/oncomatrix.R")
# source("~/maftools/R/oncoplot.R")
# pdf(paste0(mywd,"SNVs-CNVs-Pheocromocytoma.pdf"),width = 10,height = 12)
# maftools::oncoplot(maf = maf.plus.cn,genes = genes,legend_height = 10,legendFontSize = 2,annotationFontSize = 2,clinicalFeatures = c("PrimaryLocation","Tissue","SampleType"),gene_mar = 10,barcode_mar = 10)
# dev.off()

sample_annotation_data<-maf.plus.cn@clinical.data[,c("Subject.Code","PrimaryLocation","Tissue","SampleType")]
colnames(sample_annotation_data)<-c("Tumor_Sample_Barcode","PrimaryLocation","Tissue","SampleType")
sample_annotation_colors <- get_clinical_colors(sample_annotation_data)
g<-make_oncoplot(maf.plus.cn,show_sample_names = FALSE,clin_data = sample_annotation_data,clin_data_colors = sample_annotation_colors)
draw(g)
source("https://raw.githubusercontent.com/mtandon09/mt_helpers/main/helper_functions.oncoplot.R")
#source("https://raw.githubusercontent.com/mtandon09/mt_helpers/main/helper_functions.tcga.R")



make_oncoplot <- function(maf.filtered, cohort_freq_thresh = 0.01, auto_adjust_threshold=T,
                          oncomat_only=F, show_sample_names=NULL,
                          clin_data=NULL, clin_data_colors=NULL,
                          savename=NULL) {
  require(ComplexHeatmap)
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)

  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$AlteredSamples/as.numeric(maf.filtered@summary$summary[3])),
                         stringsAsFactors = F)

  ngene_max=25
  target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
  if (auto_adjust_threshold) {
    cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
  }
  ### Select genes based on the frequency threshold
  freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
  freq_genes <- freq_genes[1:ngene_max]
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
  print(genes_for_oncoplot$Hugo_Symbol)
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
      myanno <- HeatmapAnnotation(df=anno_data,col = clin_data_colors)
    }
  }

  ## Show total burden for top annotation
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  top_ha = HeatmapAnnotation("Total\nMutations" = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F),
                             annotation_name_side = "left",annotation_name_rot=90,annotation_name_gp = gpar(cex=0.7))

  # browser()

  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(cex=0.7)), show_annotation_name=F)
  # print(oncomat.plot)
  ### Make the oncoplot
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 show_pct = F,
                                 row_split=split_idx,
                                 row_title = NULL,
                                 bottom_annotation = myanno,
                                 top_annotation = top_ha,
                                 left_annotation = left_ha,
                                 show_column_names = show_sample_names,
                                 alter_fun_is_vectorized = T)#,

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

get_clinical_colors <- function(tcga_clin_data,
                                preset_columns=c("Tumor_Sample_Barcode","ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status","tissue_or_organ_of_origin")
) {
  require(maftools)
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(circlize)

  # found_columns <- intersect(colnames(tcga_clin_data), preset_columns)
  # if ( ! "Tumor_Sample_Barcode" %in% found_columns ) {
  #   stop("Clinical data must contain a 'Tumor_Sample_Barcode' column.")
  # }
  # found_columns <- setdiff(found_columns, c("Tumor_Sample_Barcode"))
  # if (length(found_columns) < 1) {
  #   stop(paste0("Clinical data must contain at least one of the columns: ", paste0(preset_columns, collapse = ", ")))
  # }

  # anno_data <- tcga_clin_data[,..found_columns]
  # anno_data$Dataset <- paste0(unique(tcga_clin_data$disease), collapse=",")

  # browser()
  found_columns<-setdiff(colnames(tcga_clin_data), c("Tumor_Sample_Barcode"))
  anno_data<-tcga_clin_data[,..found_columns]
  color_list <- lapply(found_columns, function(featurename) {

    return_val=NULL
    switch(featurename,
           "ajcc_pathologic_stage"={
             stages=sort(unique(anno_data$ajcc_pathologic_stage))
             return_val <- setNames(brewer.pal(n = length(stages), name = "Reds"), stages)
             return_val <- return_val[!is.na(names(return_val))]
           },
           "age_at_diagnosis"={
             anno_data$age_at_diagnosis <- as.numeric(as.character(anno_data$age_at_diagnosis))
             age_range=round(range(anno_data$age_at_diagnosis, na.rm = T),-1)
             age_color_length=10
             age_breaks=round(seq(age_range[1], age_range[2], length.out=age_color_length),0)
             age_color_vals=colorRampPalette(c("lightblue1","royalblue1","navy"))(age_color_length)
             return_val <- colorRamp2(age_breaks, age_color_vals)
           },
           "gender"={
             return_val <- c(female="hotpink", male="cornflowerblue")
           },
           "race"={
             races=sort(unique(anno_data$race))
             return_val <- setNames(rev(brewer.pal(n = length(races), name = "Set1")), races)
             return_val <- return_val[!is.na(names(return_val))]
           },
           "vital_status"={
             # statuses=sort(unique(anno_data$vital_status))
             return_val <- c(Alive="darkgreen",Dead="darkred")
           },
           "tissue_or_organ_of_origin"={
             tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
             return_val <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
             return_val <- return_val[!is.na(names(return_val))]
           },{
             warning(paste0("'",featurename,"' not found in columns of clinical data provided."))
             tissues=sort(unique(anno_data[[featurename]]))
             return_val <- setNames(randomcoloR::distinctColorPalette(length(tissues)), tissues)
             return_val <- return_val[!is.na(names(return_val))]

           }
    )
    return(return_val)
  })
  names(color_list) <- found_columns
  return(color_list)
}
