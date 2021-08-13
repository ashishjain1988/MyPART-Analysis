##Extract the gene named from the SF RSEM output
extractGeneNames<-function(geneNames)
{
  l<-sapply(geneNames,FUN=function(x){unlist(base::strsplit(x,split = "\\|"))})
  genesDF<-do.call(rbind,l)
  l1<-sapply(geneNames,FUN=function(x){unlist(strsplit(x,split = "\\."))})
  ensemblID<-do.call(rbind,l1)
  return(data.frame(EnsemblID=ensemblID[,1],GeneName=genesDF[,2],row.names = geneNames))
}

###Get Tissue-Specific genes from TissueEnrich
getTissueSpecificGenes<-function(rdaFilePath)
{
  #l<-load(file = paste0("~/myPART/TissueEnrichCombineExpression.rda"))
  l<-load(file = rdaFilePath)
  tsGenesPA<-dataset$`Protein-Atlas`$tissueSpecificGenes %>% filter(Group %in% c("Tissue-Enriched","Tissue-Enhanced","Group-Enriched"))
  df1 <- tsGenesPA %>% dplyr::group_by_at(vars(Gene)) %>% dplyr::summarize(Tissue = paste(Tissue,collapse = ","),Type = paste(Group,collapse = ","))
  df1<-df1[unlist(lapply(df1$Gene,FUN = function(x){return(max(dataset$`Protein-Atlas`$expressionData[x,]) >= 5)})),]
  mapping<-dataset$humanGeneMapping
  sPA<-sqldf("SELECT * FROM df1 LEFT OUTER JOIN mapping where df1.Gene == mapping.Gene")
  #genes<-dataset$humanGeneMapping %>% filter(Gene %in% tsGenes$Gene) %>% select(Gene.name)
  tsGenesGTEx<-dataset$`GTEx-Combine`$tissueSpecificGenes %>% filter(Group %in% c("Tissue-Enriched","Tissue-Enhanced","Group-Enriched"))
  df2 <- tsGenesGTEx %>% dplyr::group_by_at(vars(Gene)) %>% dplyr::summarize(Tissue = paste(Tissue,collapse = ","),Type = paste(Group,collapse = ","))
  df2<-df2[unlist(lapply(df2$Gene,FUN = function(x){return(max(dataset$`GTEx-Combine`$expressionData[x,]) >= 5)})),]
  sGTEx<-sqldf("SELECT * FROM df2 LEFT OUTER JOIN mapping where df2.Gene == mapping.Gene")
  return(list(HPA=sPA,GTEX=sGTEx))
}

###Volcano Plot for DEGs from DESeq2 Output
volcanoPlot<-function(res1,sampleNames,padj,lfc){
  #genesNames<-extractGeneNames(row.names(res1))
  #res1$GeneName<-genesNames$GeneName
  filtGenes<-res1 %>% dplyr::filter(log2FoldChange>=lfc)
  top5UpGenes<-filtGenes[order(filtGenes$padj,decreasing = F),] %>% head(n=5)
  top5UpGenes$color<-"#7570b3"
  filtGenes<-res1 %>% dplyr::filter(log2FoldChange<=-lfc)
  top5DownGenes<-filtGenes[order(filtGenes$padj,decreasing = F),] %>% head(n=5)
  top5DownGenes$color<-"#1b9e77"
  genes_label_res<-rbind(top5UpGenes,top5DownGenes)
  fontfamily="Arial"
  g<-ggplot() + geom_point(data=res1,aes(x=log2FoldChange,y=-log10(padj)),color="black") +
    geom_point(data=dplyr::filter(res1,(padj<=padj & (log2FoldChange)>=lfc)),aes(x=log2FoldChange,y=-log10(padj),color=sampleNames[1])) +
    geom_point(data=dplyr::filter(res1,(padj<=padj & (log2FoldChange)<=-lfc)),aes(x=log2FoldChange,y=-log10(padj),color=sampleNames[2])) +
    scale_color_manual("Cell Type",breaks=c(sampleNames[1],sampleNames[2]),values = c("#7570b3","#1b9e77")) +
    #theme_bw()+
    guides(color = guide_legend(override.aes = list(size=3)), size = FALSE)+
    #guides(color = TRUE, size = FALSE)+
    #scale_color(name = 'GREAT Terms',values = lables$short, labels = lables$V1) +
    labs(x=bquote(Log[2]~'(Fold Change)'), y = bquote(-Log[10]~'(P-Adjusted)')) +
    #ggtitle('Volcano plot for DEGs between ESCd and HTR-8/SVeno Cells') +
    ggrepel::geom_label_repel(data = genes_label_res,aes(x=log2FoldChange,y=-log10(padj), label = GeneName),fontface = 'bold',color = genes_label_res$color,box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50',size = 6) +
    #geom_label_repel(data = filter(res, Gene %in% mesoderm_genes),aes(x=log2FoldChange,y=-log10(padj), label = Gene),fontface = 'bold', color = 'orange',box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme_classic(base_size = 10) +
    #Axis font = 40 for ESCd VS Mesoderm
    theme(text=element_text(face = "bold"),legend.title = element_text(size=10),legend.text = element_text(size=10),axis.text = element_text(size = 20),plot.title = element_text(hjust = 0.5,size = 25),axis.title = element_text(size = 25,face = "bold"),legend.background = element_rect(colour = "black"))
  #ggsave("/Users/jain/Desktop/Figures-PDF/VolcanoPlotESCdVSMesoderm.pdf",width = 4,height = 3.5,units = "in",scale = 3,p)
  return(g)
}

###Wrapper to get the WebGestalt results
getGOEnrichmentWebGestaltR<-function(geneList,enrichDatabase,comparison,sample,parent)
{
  projectName<-paste0(comparison,"_DEGsHighIn",sample,"_",enrichDatabase)
  enrichResult<-NULL
  tryCatch({
    enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                enrichDatabase=enrichDatabase,#c("pathway_KEGG","geneontology_Biological_Process_noRedundant","pathway_Reactome"),
                                interestGene=geneList,#interestGeneFile=geneFile,
                                interestGeneType="genesymbol", referenceSet="genome", isOutput=TRUE,
                                outputDirectory=paste0(parent,"EnrichmentAnalysis"),
                                projectName=projectName,saveRawGseaResult=TRUE,
                                minNum = 5,maxNum = 2000,sigMethod = "fdr")
    projectName<-stringr::str_replace_all(projectName,pattern = "\\.",replacement = "_")
    #print(projectName)
    title<-paste0("Enriched terms from ",enrichDatabase)
    if(!is.null(enrichResult)){
      g<-getOntologyHistogram(enrichResult,title = title)
      ggsave(paste0(parent,"EnrichmentAnalysis/","Project_",projectName,"/OntologyTerms-Hist.pdf"),width = 6,height = 3,units = "in",scale = 2,g)
      g<-getOntologyDotplot(enrichResult,title = title)
      ggsave(paste0(parent,"EnrichmentAnalysis/","Project_",projectName,"/OntologyTerms-Dotplot.pdf"),width = 6,height = 3,units = "in",scale = 2,g)
      if(enrichDatabase == "geneontology_Biological_Process_noRedundant")
      {
        saveGOTermsTreemap(enrichResult,paste0(parent,"EnrichmentAnalysis/","Project_",projectName,"/GOTerms-treemap.pdf"))
      }
    }},
    error=function(error)
    {
      message("Original error message:")
      message(error)
    },finally = function()
    {
      return(enrichResult)
    })
}

###Histogram to show the top 10 GO results from the WEbGestalt output
library(stringr)
getOntologyHistogram<-function(enrichResult,pValueCutOff=0.05,foldChangeCutoff=2,nTerms=10,title=""){
  enrichResultTop<-enrichResult[order(enrichResult$FDR),]
  if(nrow(enrichResult) >=nTerms)
  {
    enrichResultTop<-enrichResultTop[1:nTerms,]
  }
  enrichResultTop <- enrichResultTop[,c(2,7,9)]
  enrichResultTop$logFDR <- -log10(enrichResultTop$FDR)
  p<-ggplot(enrichResultTop,aes(x=reorder(description,-c(1:nrow(enrichResultTop))),y=logFDR))+
    geom_bar(stat = 'identity', position = position_dodge(width=1),fill="black")+
    coord_flip()+
    theme_bw()+
    guides(fill = FALSE)+
    #scale_y_continuous(breaks = 1:8, labels = c(1:3,"Break",7,"Break",12:13)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme(text=element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20),axis.text = element_text(size=12),axis.title = element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(x='', y = bquote(-Log[10]~'(FDR)'))+
    ggtitle(title)
  return(p)
}

###Dotplot to show the top 10 GO results from the WEbGestalt output
getOntologyDotplot<-function(enrichResult,pValueCutOff=0.05,foldChangeCutoff=2,nTerms=10,title=""){
  enrichResultTop<-enrichResult[order(enrichResult$FDR),]
  if(nrow(enrichResult) >=nTerms)
  {
    enrichResultTop<-enrichResultTop[1:nTerms,]
  }
  enrichResultTop <- enrichResultTop[,c(2,5,7,9)]
  enrichResultTop$logFDR <- -log10(enrichResultTop$FDR)
  enrichResultTop$GeneOverlap<-enrichResultTop$overlap

  p <- ggplot(enrichResultTop, aes(x=logFDR, y=reorder(description,-c(1:nrow(enrichResultTop))), size=GeneOverlap, color=enrichmentRatio)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "Enrichment ratio") + #,guide=guide_colorbar(reverse=TRUE)
    #scale_y_discrete(labels = label_func) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme(text=element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20),axis.text = element_text(size=12),axis.title = element_text(size=15)) +
    labs(y='', x = bquote(-Log[10]~'(FDR)'))+
    ggtitle(title)
  #theme_dose(font.size) +
  #scale_size(range=c(3, 8))
  return(p)
}

library(rrvgo)
saveGOTermsTreemap<-function(enrichResult,fileName)
{
  simMatrix <- calculateSimMatrix(enrichResult$geneSet,orgdb="org.Hs.eg.db",ont="BP",method="Rel")
  minValue<-ifelse(min(enrichResult$FDR[enrichResult$FDR !=0]) < 10e-7,min(enrichResult$FDR[enrichResult$FDR !=0]),10e-7)
  scores <- setNames(-log10(enrichResult$FDR+(minValue)), enrichResult$geneSet)
  reducedTerms <- reduceSimMatrix(simMatrix,scores,threshold=0.7,orgdb="org.Hs.eg.db")
  pdf(fileName)
  treemapPlot(reducedTerms)
  dev.off()
}

getGeneEffectsInDepMap<-function(genes=NULL,depMapVersion=NULL,groupBySubtypes=FALSE,file_path=NULL){

  if (is.null(genes)) {
    stop("Need to enter the gene or gene list")
  }
  require(depmap)
  require("ExperimentHub")
  require(tidyverse)
  require(gridExtra)
  require(cowplot)
  eh <- ExperimentHub()
  dep<-query(eh, "depmap")
  ##Get Crispr dependency Scores
  ##Version 21Q1
  pedDepMapDataInfo<-read.table("~/myPART/DepMapPedSampleInfo.txt",sep = "\t",header = T)
  pedDepMapDataInfo <- pedDepMapDataInfo[pedDepMapDataInfo$Class_for_Manuscript == "Pediatric",]
  crisprData<-eh[["EH5358"]]#depmap::depmap_crispr()
  metaData<-eh[["EH5362"]]#depmap::depmap_metadata()
  t<-data.frame(crisprData %>% dplyr::filter(depmap_id %in% pedDepMapDataInfo$DepMap_ID))
  plotsLists<-list()
  for(gene in genes)
  {
    #gene<-"CHEK1"
    #d$gene<-unlist(lapply(d$gene,FUN=function(x){return(unlist(str_split(x," "))[1])}))
    d <- t[t$gene_name %in% gene,]
    d1 <- merge(d,metaData,by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name,primary_disease,subtype_disease,dependency,cell_line.x)
    d1$subtype_disease[is.na(d1$subtype_disease)]<-d1$primary_disease[is.na(d1$subtype_disease)]
    d1$subtype_disease[d1$subtype_disease == "Undifferentiated"]<-paste0(d1$primary_disease[d1$subtype_disease == "Undifferentiated"],"(",d1$subtype_disease[d1$subtype_disease == "Undifferentiated"],")")
    medianDep<-d1 %>% group_by(subtype_disease) %>% summarise(Mean=mean(dependency),Median=median(dependency))
    #print(paste0(genes,"-",round(min(medianDep$Median),digits = 2)))
    g<-ggplot2::ggplot(d1,aes(x=subtype_disease,y=dependency))+
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 10),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
      ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
      ggplot2::labs(x="", y = "Normalized Dependency Score") +
      ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
      #geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5) +
      coord_flip()
    g1<-ggplot2::ggplot(d1,aes(x=dependency)) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 12),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
      ggplot2::labs(x="", y = "") +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept=-1, linetype="dashed", color = "red")
    p <- cowplot::plot_grid(g1, g, align = "v",ncol = 1,rel_heights=c(1,2))
    #save_plot(paste0(file_path,"Demap-PedCellLines-",round(min(medianDep$Median),digits = 2),"-CRISPR-",gene,".pdf"), p,base_height = 10,base_width = 8)
    plotsLists[[gene]]<- list(ggHistogram=g1,ggBoxplot=g,combined=p)
  }
  return(plotsLists)
}

generateVariantTableFromMAFData <- function(maf, use_syn=FALSE, extra_cols=c()) {

  output_data<-maf
  if (all(c("Tumor_Seq_Allele1","Tumor_Seq_Allele2") %in% colnames(output_data))) {
    output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
  }
  if (all(c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2") %in% colnames(output_data))) {
    output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")
  }

  if (! "tumor_freq" %in% colnames(output_data)) {
    # browser()
    if (all(c("t_depth","t_alt_count")%in% colnames(output_data))) {
      output_data$tumor_freq <- as.numeric(as.character(output_data$t_alt_count))/as.numeric(as.character(output_data$t_depth))
    }
  }
  cols_for_table <- c("Hugo Symbol" = "Hugo_Symbol",
                      "Sample ID" = "Tumor_Sample_Barcode",
                      "Variant Classification"="Variant_Classification",
                      "Variant Type"="Variant_Type",
                      "Consequence"="Consequence",
                      "Chromosome"="Chromosome","Start Position" ="Start_Position","End Position"="End_Position","Strand"="Strand",
                      "Reference Allele"="Reference_Allele",
                      "Tumor Genotype"="tumor_genotype",
                      "Normal Genotype"="normal_genotype",
                      "Known Effects ClinVar"="CLIN_SIG",
                      "Transcript Change"="HGVSc",
                      "Protein Change"="HGVSp_Short",
                      "Normal Depth"="n_depth",
                      "Normal Ref Depth"="n_ref_count",
                      "Normal Alt Depth"="n_alt_count",
                      "Tumor Depth"="t_depth",
                      "Tumor Ref Depth"="t_ref_count",
                      "Tumor Alt Depth"="t_alt_count",
                      "Tumor Alt Frequency"="tumor_freq",
                      "Existing Annotation"="Existing_variation",
                      "gnomAD Frequency"="gnomAD_AF",
                      "ExAC Frequency"="ExAC_AF",
                      "1000Genomes Frequency"="AF",
                      "Effect Prediction - SIFT"="SIFT",
                      "Effect Prediction - PolyPhen"="PolyPhen"
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


