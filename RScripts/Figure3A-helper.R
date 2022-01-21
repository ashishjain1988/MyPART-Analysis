diagnosisColorRNASeq <- c("#78DD9E","#A54EE1","#D2B4C8","#D7815D","#87D2DC","#BCE25A","#D5DAA9","#DC6DBF","#8688D5")
names(diagnosisColorRNASeq)<-c("Adrenocortical carcinoma","Anaplastic thyroid carcinoma","Medullary thyroid carcinoma","Neuroendocrine-small intestine","Neuroendocrine-pancreas","Neuroendocrine-lung","Papillary thyroid carcinoma","Pheochromocytoma","Normal")

diagnosisColorDNASeq <- c("#78DD9E","#A54EE1","#D2B4C8","#D7815D","#D5DAA9","#DC6DBF","#8688D5")
names(diagnosisColorDNASeq)<-c("Adrenocortical carcinoma","Anaplastic thyroid carcinoma","Medullary thyroid carcinoma","Neuroendocrine","Papillary thyroid carcinoma","Pheochromocytoma","Normal")

tissueColor<-c("#DC5ECD","#806AD9","#D2DC98","#DA6074","#D49555","#CDE5D5","#B5C1DF","#D3ABAB","#73D9AC","#9E3DE7","#69D1DE","#D298D3","#688DD0","#85E66D","#D4DE54")
names(tissueColor)<-c("Adrenal" ,"Colon","Kidney","Liver","Lung","Lymph Node","Neck","Ovary","Pancreas","Pelvis","Scalp","Small Intestine","Spinal Cord","Stomach","Thyroid")

tissueColorDNASeq<-c("#DC5ECD","#806AD9","#D2DC98","#DA6074","#D49555","#CDE5D5","#B5C1DF","#D3ABAB","#73D9AC","#9E3DE7","#69D1DE","#D298D3","#688DD0","#85E66D","#D4DE54","#799D74","#9E3DE7")
names(tissueColorDNASeq)<-c("Adrenal" ,"Colon","Kidney","Liver","Lung","Lymph Node","Neck","Ovary","Pancreas","Pelvis","Scalp","Intestine","Spinal Cord","Stomach","Thyroid","Rectum","Ileum")

tumorSiteColor<-c("Metastasis"="black","Primary Site"="#BB70CB","Recurrence"="#B8D97B")
genderColor<-c("Male"="turquoise","Female"="brown")
raceColor<-c("White"="#B35806","Unknown"="#FDBC6B","Other"="black","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue")
getAgeColor<-function(ageVector){colorRamp2(c(0,max(ageVector,na.rm = T)), c("white","blue"))}
getTMB500Color<-function(TMBValues){colorRamp2(c(0,max(TMBValues,na.rm = T)), c("white","red"))}

diagnosisShortNames <- c("ACC","ATC","MTC","NET","NET-SI","NET-pancreas","NET-lung","PTC","PCC","Normal")
names(diagnosisShortNames)<-c("Adrenocortical carcinoma","Anaplastic thyroid carcinoma","Medullary thyroid carcinoma","Neuroendocrine","Neuroendocrine-small intestine","Neuroendocrine-pancreas","Neuroendocrine-lung","Papillary thyroid carcinoma","Pheochromocytoma","Normal")

diagnosisColorRNASeqShort <- c("#78DD9E","#A54EE1","#D2B4C8","#D7815D","#87D2DC","#BCE25A","#D5DAA9","#DC6DBF","#8688D5","#D7815D")
names(diagnosisColorRNASeqShort)<-c("ACC","ATC","MTC","NET-SI","NET-pancreas","NET-lung","PTC","PCC","Normal","NET")

### Creates matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = unique(g))]
  }

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                #x = unique(as.character(x))
                                x = as.character(x)
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                  # xvc = paste0(xvc, collapse="|")
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}


### Define the colors for different annotations
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
           },"Age"={
             # tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
             # return_val <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
             return_val <- colorRamp2(c(0,max(anno_data$Age,na.rm = T)), c("white","blue"))
           },"TMB_TSO500"={
             # tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
             # return_val <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
             return_val <- colorRamp2(c(0,max(anno_data[[featurename]],na.rm = T)), c("white","red"))
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


### Define colors for mutation types
mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
                     In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
                     In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
                     Amp="green2",Del="darkred",
                     no_variants="#d6d6d6", Pathogenic="black",VUS="grey50")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))

### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
  },
  "Amp" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
              gp = gpar(fill = mutation_colors["Amp"], col = NA))
  },
  "Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
              gp = gpar(fill = mutation_colors["Del"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "Pathogenic" = function(x, y, w, h) {
    # grid.points(x, y, pch = 18, size=w, gp=gpar(col=col["pathogenic"]))
    # grid.rect(x, y, w*0.7, h*0.2,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    # grid.rect(x, y, w*0.1, h*0.7,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["Pathogenic"], fill = NA, lwd=5))
  },
  "VUS" = function(x, y, w, h) {
    # grid.points(x, y, pch = 3, size=w,gp=gpar(col=col["VUS"], lwd=3))
    # grid.rect(x, y, w*0.2, h-unit(0.5, "mm"),
    #           gp = gpar(fill = col["VUS"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["VUS"], fill = NA, lwd=5))
  }
)


