library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)
library(ggplot2)
library(dplyr)

###########Figure 1A###################
source("~/myPART/MyPART-Analysis/RScripts/helperFunctions.R")
source("~/myPART/MyPART-Analysis/RScripts/Figure3A-helper.R")
parent<-"~/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/"
l1<-load(file = paste0(parent,"vsdNormalizedCounts.rda"))
features<-featuresFilt
vsdNormalizedCounts <- vsdNormalizedCounts[,row.names(features)]
features$Tissue <- tools::toTitleCase(features$Tissue)
features$TumorSite <- tools::toTitleCase(features$TumorSite)
pca <- prcomp(t(as.matrix(vsdNormalizedCounts)),scale. = T)
eigs <- pca$sdev^2
var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%",sep = "")
pcaData<-pca$x
colnames(pcaData)<-paste(colnames(pcaData),rep("(",ncol(pcaData)),var_explained,rep(")",ncol(pcaData)),sep = "")

pca_exp<-pcaData[,1:5]
featuresFilt<-features#features[grepl("Normal|Gastrointestinal",features$New.Diagnosis),]
hm_data<-as.matrix(t(pca_exp[row.names(featuresFilt),]))
column_annotations = HeatmapAnnotation(df = featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = getAgeColor(featuresFilt$Age)#colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","blue"))
#col_rin = colorRamp2(c(0,max(featuresFilt$RIN)), c("white","yellow"))
#n <- distinctColorPalette(length(levels(as.factor(featuresFilt$Tissue))))
#names(n)<-(levels(as.factor(featuresFilt$Tissue)))
col_tissue <-tissueColor[featuresFilt$Tissue]

#n1 <- distinctColorPalette(length(levels(as.factor(featuresFilt$New.Diagnosis))))
#names(n1)<-(levels(as.factor(featuresFilt$New.Diagnosis)))
#n1 <- c("#78DD9E","#A54EE1","#D2B4C8","#D7815D","#87D2DC","#BCE25A","#D5DAA9","#DC6DBF","#8688D5")
#names(n1)<-unique(featuresFilt$New.Diagnosis)
col_disease <-diagnosisColorRNASeq[featuresFilt$New.Diagnosis]

df<-featuresFilt[,c("New.Diagnosis","Tissue","Age","Sex","Race","TumorSite")]
colnames(df)<-c("Diagnosis","Tissue","Age","Sex","Race","TumorSite")
#col_tissue <- distinctColorPalette()
ha = HeatmapAnnotation(
  #df = features[,c(5,6,7,13,15)],
  #df = featuresFilt[,c(3,4,5,6,13,7)],
  df = df,
  col = list(
             Diagnosis = col_disease,
             Tissue = col_tissue,#c("Normal" = "black","Adrenocortical carcinoma" = "red"),#col_cluster,#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             Age = col_age,
             Sex = genderColor,#c("Male"="turquoise","Female"="brown"),
             #RIN = col_rin,
             Race =raceColor,#c("White"="#B35806","Unknown"="#FDBC6B","Other"="black","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue"),
             TumorSite=tumorSiteColor#c("Metastasis"="turquoise","Primary Site"="brown","Recurrence"="black")
             #cluster = col_cluster
  ),
  #annotation_height = unit(rep(6,6), rep("mm",6)),
  #height = unit(40,"mm"),
  annotation_name_gp =  gpar(fontsize = 14,fontface = 2),
  #annotation_legend_param = list(title_gp = gpar(fontsize = 15, fontface = 2),labels_gp = gpar(fontsize = 12)),
  #annotation_height = unit(rep(2,6), "cm"),
  height=unit(5*0.4,'inches'),
  annotation_legend_param=list(labels_gp=gpar(fontsize = 12),
                               #grid_height=unit(3,"mm"),
                               #grid_width=unit(3,"mm"),
                               # ncol = 5,
                               nrow = 4,
                               legend_direction = "vertical",
                               title_gp=gpar(fontsize=15,
                                             # legend_height=unit(0.3,"npc"),
                                             fontface="bold"))
  #annotation_height = c(TumorSite = unit(8, "mm"),Diagnosis=unit(8, "mm"),Age=unit(8, "mm"),Sex=unit(8, "mm"),Race=unit(8, "mm"),TumorSite=unit(8, "mm")), height = unit(48, "mm")
)
# ha@height<-unit(48,"mm")
# ha@anno_size<-unit(rep(8,6), "mm")

#ha1<-re_size(ha,height = 40)
cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(pca_exp), na.rm = T), max(abs(pca_exp), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=FALSE,
                                    #column_names_rot = 45,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    clustering_distance_columns = "euclidean",
                                    clustering_method_columns = "complete",
                                    row_names_gp = gpar(fontsize = 14,fontface = 2),
                                    #column_names_gp = gpar(fontsize = 14,fontface = 2),
                                    show_heatmap_legend = F,
                                    heatmap_height = unit(100,"mm")
                                    ) # Turning off to control the placement

pdf(paste0(parent,"PCsSampleClustering-Subtypes-All.pdf"),width = 7,height = 8)
ComplexHeatmap::draw(cheatmap, show_annotation_legend = TRUE,annotation_legend_side = "bottom")
dev.off()

###########Figure 1B###################
f<-data.frame(Diagnosis=as.numeric(as.factor(features$New.Diagnosis)))
f$Sex<-as.numeric(as.factor(features$Sex))
f$Race<-as.numeric(as.factor(features$Race))
f$Age<-(features$Age)
cormatrix<-round(cor(pcaData[,1:5],f, method = "spearman",use = "na.or.complete"),2)
cheatmap <- ComplexHeatmap::Heatmap(t(cormatrix),
                                    col=colorRamp2(seq(-1, 1, length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                      grid.text(t(cormatrix)[i, j], x, y,gp=gpar(fontsize = 20,fontface = 2))
                                    },
                                    bottom_annotation = NULL,
                                    show_column_names=T,
                                    column_names_rot = 45,
                                    cluster_rows = T,
                                    cluster_columns = F,
                                    row_names_gp = gpar(fontsize = 20,fontface = 2),
                                    column_names_gp = gpar(fontsize = 20,fontface = 2),
                                    show_heatmap_legend = FALSE)
pdf(paste0(parent,"PCsFeaturesCorrelation.pdf"),width = 10,height = 8)
draw(cheatmap, show_annotation_legend = TRUE)
dev.off()

###########Figure 1C###################
#parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/"
o<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/","cancerSpecificGenes-TPM.rda"))
drugTargets<-read.table(paste0("/Users/jaina13/myPART/AllSamplesPipeliner/","CancerSpecificGenes/DrugBankTargets-Format.tsv"),header = TRUE,sep = "\t")
drugBankGeneMapping<-read.table(paste0("/Users/jaina13/myPART/AllSamplesPipeliner/","CancerSpecificGenes/genesDrugBank.tsv"),header = TRUE,sep = "\t")
drugTargetsMap<-merge(x=drugTargets,y=drugBankGeneMapping,by.x="targetID",by.y="gene_claim_name",all.x=TRUE,all.y=FALSE)
#teGenes<-getTissueSpecificGenes(rdaFilePath = "~/myPART/TissueEnrichCombineExpression.rda")
proteinCodingGenes<-read.table(paste0("~/myPART/FusionAnalysis/hg38.proteinCodingGenes.Ensembl.txt"),sep = "\t",header = T)
l<-list()
for(tissue in unique(features$New.Diagnosis))
{
  if(!grepl("Normal",tissue))
  {
    c<-cellSpecificGenesOutput %>% dplyr::filter(Group %in% c("Tissue-Enriched")) %>% dplyr::filter(Tissue == tissue)
    c2<-sapply(c$Gene,FUN = function(x){
      return(paste0(drugTargetsMap[drugTargetsMap$gene_name %in% x,2],collapse = ","))
    })
    c$drugBank<-c2
    c<-c[c$Gene %in% proteinCodingGenes$Gene.name,]
    l[[str_replace_all(tissue,pattern = "[ /]",replacement = ".")]]<-c
  }
}
#o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
csGenes<-do.call(rbind,l)
csGenes$isDrug<-"No"
csGenes[csGenes$drugBank != "",]$isDrug<-"Yes"
#csGenes<-csGenes[csGenes$isDrug =="Yes",]
g<-ggplot(csGenes, aes(fill=isDrug, x=(Tissue))) +
  #g<-ggplot(csGenes, aes(x=(Tissue))) +
  geom_bar(position="stack") +
  scale_fill_manual(values = c("grey","blue")) +
  guides(fill=guide_legend(title="isDrugTarget")) +
  labs(x="", y = "# Cancer-Specific Genes") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.title = element_text(size=14,face = "bold"),legend.text = element_text(size=14),axis.title = element_text(size = 14,face = "bold"),axis.text = element_text(size = 16,face = "bold",colour = "black")) +
  theme(plot.margin = unit(c(1,1,1,2), "cm")) +
  #margin(t, r, l, b)
  xlab("")
ggsave(paste0(parent,"EndocrineSubgroupResults/cancerSpecificGenesBarplot-DrugTarget-AllSamples.pdf"),plot = g,height = 8,width = 10)

###########Figure 1D [cancer-specific genes HeatMap]###################
csGenesDrugBank<-csGenes[csGenes$isDrug =="Yes",]
tpmValuescsGenes<-tpmValuesFull1[csGenesDrugBank$Gene,]
tpmValuescsGenes <- log2(tpmValuescsGenes+1)
# r<-data.frame(t(apply(tpmValuescsGenes, 1, scale)))
# colnames(r)<-colnames(tpmValuescsGenes)
# tpmValuescsGenes<-r

features$Tissue <- tools::toTitleCase(features$Tissue)
features$TumorSite <- tools::toTitleCase(features$TumorSite)
featuresFilt<-features#features[grepl("Normal|Gastrointestinal",features$New.Diagnosis),]
featuresFilt[grepl("Normal",featuresFilt$New.Diagnosis),]$New.Diagnosis<-"Normal"
hm_data<-as.matrix((tpmValuescsGenes[,row.names(featuresFilt)]))
column_annotations = HeatmapAnnotation(df = featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = getAgeColor(featuresFilt$Age)#colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","blue"))
col_tissue <-tissueColor[featuresFilt$Tissue]
col_disease <-diagnosisColorRNASeq[featuresFilt$New.Diagnosis]

df<-featuresFilt[,c("New.Diagnosis","Tissue","Age","Sex","Race","TumorSite")]
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
  #annotation_height = unit(rep(2,6), "cm")
  #simple_anno_size = unit(1, "cm"), height = unit(6, "cm")
)

#ha1<-re_size(ha,height = 15,annotation_height = unit(8, "mm"))
#ha_row = rowAnnotation(foo = anno_mark(at = c(1,3:8,26,28,31,45:46), labels = row.names(tpmValuescsGenes)[c(1,3:8,26,28,31,45:46)]))
#subgroup = sample(letters[1:3], 100, replace = TRUE, prob = c(1, 5, 10))
library(ggplotify)
rg = range(hm_data)
panel_fun = function(index, nm) {
  pushViewport(viewport(xscale = rg, yscale = c(0, 2)))
  grid.rect()
  #grid.xaxis(gp = gpar(fontsize = 8))
  #grid.boxplot(hm_data[index, ], pos = 1, direction = "horizontal")
  #png_obj <- image_read(paste0(parent,"TPH1-expBoxplot.pdf"))
  #png_obj <- image_trim(png_obj)
  #grid.raster(png_obj,width=unit(1,"npc"), height=unit(1,"npc"))
  #print(index)
  gene <- row.names(hm_data[index, ,drop=FALSE])
  d<-reshape2::melt(tpmValuesFull1[gene, ])
  d$variable <- colnames(tpmValuescsGenes)
  f<-featuresFilt[,c("RTNo","PrimaryCancerTissue","New.Diagnosis","TumorSite")]
  #f$New.Diagnosis[f$New.Diagnosis == "Adrenocortical carcinoma-recurrence"]<-"Adrenocortical carcinoma-metastasis"
  f$samples<-row.names(f)
  final_data <- sqldf::sqldf("SELECT * FROM d LEFT OUTER JOIN f where d.variable == samples")
  #print(final_data)
  g<-ggplot(final_data,aes(x=New.Diagnosis,y=value,fill = New.Diagnosis,color = New.Diagnosis))+
    theme_classic(base_size = 10) +
    scale_fill_manual(breaks=unique(names(col_disease)),values=diagnosisColorRNASeq[unique(names(col_disease))])+
    scale_color_manual(breaks=unique(names(col_disease)),values=diagnosisColorRNASeq[unique(names(col_disease))])+
    #theme(text=element_text(face = "bold",colour = "black"),title = element_text(size = 16),axis.text = element_text(size = 20,colour = "black"),axis.title = element_text(size = 20,face = "bold"),legend.background = element_rect(colour = "black")) +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.2)) +
    theme(plot.title = element_text(hjust = 0.5,size = 8,face = "bold"),axis.text = element_text(size = 4,colour = "black"),axis.title = element_text(size = 5,colour = "black"))+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1)+
    #geom_point(size=1.5)+
    guides(fill="none",color="none") +
    labs(x="", y = "TPM") +
    ggtitle(paste0(gene))
    #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
    #theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
    #coord_flip()
  grid.draw(as.grob(g))
  popViewport()
}
genes<-which(row.names(hm_data) %in% c("GNRHR","TPH1","NEUROD1","CALCB","MMP9","SLC6A2","DRD2"))
#names(genes)<-paste0("Gene",1:6)
anno = anno_zoom(align_to = as.list(genes), which = "row", panel_fun = panel_fun,
                 size = unit(2, "cm"), gap = unit(1, "cm"), width = unit(4, "cm"))

cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(tpmValuescsGenes), na.rm = T), max(abs(tpmValuescsGenes), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=FALSE,
                                    #column_names_rot = 45,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    clustering_distance_columns = "euclidean",
                                    clustering_method_columns = "complete",

                                    row_names_side = "left",
                                    show_row_names = FALSE,
                                    right_annotation = rowAnnotation(foo1=anno),
                                    row_names_gp = gpar(fontsize = 10,fontface = 2),
                                    #left_annotation =ha_row,
                                    left_annotation = rowAnnotation(Reason = csGenesDrugBank$Tissue, col=list(Reason=diagnosisColorRNASeq[unique(csGenesDrugBank$Tissue)]), annotation_width = unit(0.3, "mm"), show_annotation_name=F,show_legend=FALSE),

                                    #column_names_gp = gpar(fontsize = 14,fontface = 2),
                                    show_heatmap_legend = F
) # Turning off to control the placement

pdf(paste0(parent,"PCsSampleClustering-csGenes.pdf"),width = 12,height = 7)
ComplexHeatmap::draw(cheatmap, show_annotation_legend = TRUE)
dev.off()


##############GO Term barplot for cancer-specific genes##########
#parent<-""
modulesGOTerms<-c("Pheochromocytoma","Papillary_thyroid_carcinoma","Anaplastic_thyroid_carcinoma")
names(modulesGOTerms)<-c("Pheochromocytoma","Papillary thyroid carcinoma","Anaplastic thyroid carcinoma")
#####Read top 5 GO terms of each of the four modules
top5TermsAll<-list()
maxTerms<-5
for(co in modulesGOTerms)
{
  data<-read.table(paste0(parent,"CancerSpecificGenes/EnrichmentAnalysis/Project_",co,"_geneontology_Biological_Process_noRedundant/enrichment_results_",co,"_geneontology_Biological_Process_noRedundant.txt"),sep = "\t",header = T)
  enrichmentRes<-data[,c(1,2,7,9)] %>% dplyr::filter(FDR <= 0.05)
  top5terms<-enrichmentRes[order(enrichmentRes$enrichmentRatio,decreasing = TRUE),][1:min(nrow(enrichmentRes),maxTerms),]
  top5terms$Diagnosis<-names(modulesGOTerms[modulesGOTerms == co])
  top5TermsAll[[co]]<-top5terms
}
top5TermsAllDF<-do.call("rbind",top5TermsAll)
top5TermsAllDF$description[6]<-paste0(top5TermsAllDF$description[6]," ")
top5TermsAllDF$description <- tools::toTitleCase(top5TermsAllDF$description)
#data[6,"V3"]<-data[6,"V3"]-5
#data[7,"V3"]<-data[7,"V3"]-2
#top5TermsAllDF$Module<-as.factor(top5TermsAllDF$Module,levels=c("M4", "M5", "M15","M16"))
p<-ggplot(top5TermsAllDF,aes(x=reorder(description,-c(1:nrow(top5TermsAllDF))),y=enrichmentRatio,fill=Diagnosis))+
  geom_bar(stat = 'identity', position = position_dodge(width=1))+
  #scale_fill_manual(labels = c("M4", "M5", "M15","M16"), values = )+
  #scale_fill_discrete(breaks=levels(top5TermsAllDF$Module)) +
  scale_fill_manual(breaks=c("Anaplastic thyroid carcinoma", "Papillary thyroid carcinoma","Pheochromocytoma"),
                    values=c("Anaplastic thyroid carcinoma"="#A54EE1","Papillary thyroid carcinoma"="#D5DAA9","Pheochromocytoma"="#DC6DBF"))+
  coord_flip()+
  theme_bw()+
  # geom_hline(yintercept = 3.9,linetype="solid") +
  # geom_hline(yintercept = 4.1,linetype="solid") +
  # geom_hline(yintercept = 5.9,linetype="solid") +
  # geom_hline(yintercept = 6.1,linetype="solid") +
  #guides(fill = "none")+
  #scale_y_continuous(breaks = 1:8, labels = c(1:3,"Break",7,"Break",12:13)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
  theme(text=element_text(size = 15),legend.title = element_text(size=18,colour = "black",face = "bold"),legend.text = element_text(size=18,colour = "black"),plot.title = element_text(hjust = 0.5,size = 20),axis.text = element_text(size=18,colour = "black"),axis.title = element_text(size=20,colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #labs(x='', y = bquote(-Log[10]~'(FDR)'))
  labs(x='', y = "Enrichment Ratio")
#facet_wrap(~Module,ncol = 4)
#ggtitle('Module Top 5 GO Biological Process terms')
ggsave(paste0(parent,"csGenes.GOterms.pdf"),width = 6.6,height = 3.7,units = "in",scale = 2,p)


###########Figure 1D [Expression Boxplots]###################
o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
#genes<-c("GNRHR","MMP","TPH1","CELA2A","KCNJ6","F5","GHRHR")
gene<-"NEUROD1"
tpmvaluesGene<-tpmValuesFull1[gene,]
d<-reshape2::melt(tpmvaluesGene)
f<-features[,c("RTNo","PrimaryCancerTissue","New.Diagnosis","TumorSite")]
#f$New.Diagnosis[f$New.Diagnosis == "Adrenocortical carcinoma-recurrence"]<-"Adrenocortical carcinoma-metastasis"
f$samples<-row.names(f)
final_data <- sqldf::sqldf("SELECT * FROM d LEFT OUTER JOIN f where d.variable == samples")
diagnosisColorRNASeq1 <- c("#78DD9E","#A54EE1","#D2B4C8","#D7815D","#87D2DC","#BCE25A","#D5DAA9","#DC6DBF",rep("#8688D5",8))
names(diagnosisColorRNASeq1)<-c("Adrenocortical carcinoma","Anaplastic thyroid carcinoma","Medullary thyroid carcinoma","Neuroendocrine-small intestine","Neuroendocrine-pancreas","Neuroendocrine-lung","Papillary thyroid carcinoma","Pheochromocytoma","Normal-adrenal","Normal-kidney","Normal-lung","Normal-pancreas","Normal-small intestine","Normal-spinal cord","Normal-stomach","Normal-thyroid")
col_disease1<-diagnosisColorRNASeq1[f$New.Diagnosis]

g<-ggplot(final_data,aes(x=New.Diagnosis,y=value,fill = New.Diagnosis,color = New.Diagnosis))+
  theme_classic(base_size = 10) +
  theme(text=element_text(face = "bold",colour = "black"),title = element_text(size = 16),axis.text = element_text(size = 20,colour = "black"),axis.title = element_text(size = 20,face = "bold"),legend.background = element_rect(colour = "black")) +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.2)) +
  scale_fill_manual(breaks=unique(names(col_disease1)),values=diagnosisColorRNASeq1[unique(names(col_disease))])+
  scale_color_manual(breaks=unique(names(col_disease1)),values=diagnosisColorRNASeq1[unique(names(col_disease))])+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE)+
  #geom_point(size=1.5)+
  guides(size="none",fill="none",color="none") +
  labs(x="", y = "Transcripts Per Million (TPM)") +
  ggtitle(paste0("Expression of ",gene, "\n(Neuronal Differentiation Factor 1)"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  coord_flip()
ggsave(paste0(parent,"/",gene,"-expBoxplot.pdf"),plot = g,height = 10,width = 12)

###########Figure 1D [DepMap Boxplots]###################
require(depmap)
require("ExperimentHub")
require(tidyverse)
require(gridExtra)
require(cowplot)
eh <- ExperimentHub()
dep<-query(eh, "depmap")
##Get Crispr dependency Scores
##Version 21Q1
#pedDepMapDataInfo<-read.table("~/myPART/DepMapPedSampleInfo.txt",sep = "\t",header = T)
#pedDepMapDataInfo <- pedDepMapDataInfo[pedDepMapDataInfo$Class_for_Manuscript == "Pediatric",]
crisprData<-eh[["EH5358"]]#depmap::depmap_crispr()
metaData<-eh[["EH5362"]]#depmap::depmap_metadata()

#genes<-hubGenes$gene[1:10]#c("STAR","CYP17A1","CYP21A2","FDXR","CYP11A1","POR","CYP11B1","NR5A1","FDX1","NR0B1")#c("AURKB","SMARCA4")
#genes<-c("FGF4","PLAT","CCKBR","KCNA1","KIT","GRIN3A")
#genes<-c("PLXNA2","TNRC6C" ,"LTBP1" , "PJA2" ,  "PTCH1"  ,"USO1")
gene<-"TPH1"
#t<-data.frame(crisprData %>% dplyr::filter(depmap_id %in% pedDepMapDataInfo$DepMap_ID))
t<-data.frame(crisprData)
plotsLists<-list()
#for(gene in genes)
{
  #gene<-"CHEK1"
  #d$gene<-unlist(lapply(d$gene,FUN=function(x){return(unlist(str_split(x," "))[1])}))
  d <- t[t$gene_name %in% gene,]
  d1 <- merge(d,metaData,by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name,primary_disease,subtype_disease,dependency,cell_line.x)
  d1$subtype_disease[is.na(d1$subtype_disease)]<-d1$primary_disease[is.na(d1$subtype_disease)]
  d1$subtype_disease[d1$subtype_disease == "Undifferentiated"]<-paste0(d1$primary_disease[d1$subtype_disease == "Undifferentiated"],"(",d1$subtype_disease[d1$subtype_disease == "Undifferentiated"],")")
  medianDep<-d1 %>% group_by(subtype_disease) %>% summarise(Mean=mean(dependency),Median=median(dependency))
  #print(paste0(genes,"-",round(min(medianDep$Median),digits = 2)))
  g<-ggplot2::ggplot(d1,aes(x=primary_disease,y=dependency))+
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 20,colour = "black"),axis.title = element_text(size = 20,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
    #ggplot2::geom_point(size=1.5)+
    #ggplot2::guides(size="none") +
    ggplot2::labs(x="", y = "Normalized Dependency Score") +
    ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
    theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
    #geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5) +
    coord_flip()
  # g1<-ggplot2::ggplot(d1,aes(x=dependency)) +
  #   ggplot2::theme_classic(base_size = 10) +
  #   ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 12),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
  #   ggplot2::labs(x="", y = "") +
  #   ggplot2::geom_histogram() +
  #   ggplot2::geom_vline(xintercept=-1, linetype="dashed", color = "red")
  # p <- cowplot::plot_grid(g1, g, align = "v",ncol = 1,rel_heights=c(1,2))
  #save_plot(paste0(parent,"Demap-PedCellLines-",round(min(medianDep$Median),digits = 2),"-CRISPR-",gene,".pdf"), p,base_height = 10,base_width = 8)
  #plotsLists[[gene]]<- list(ggHistogram=g1,ggBoxplot=g,combined=p)
  plotsLists[[gene]]<- list(ggBoxplot=g)
}
ggsave(paste0(parent,"EndocrineSubgroupResults/",gene,"-DepMapBoxplot.pdf"),plot = plotsLists[[gene]]$ggBoxplot,height = 10,width = 12)

