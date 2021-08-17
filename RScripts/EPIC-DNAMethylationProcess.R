BiocManager::install("minfi")
BiocManager::install("ChAMP")
library(minfi)
library(ChAMP)
library(openxlsx)
baseDir <- "/Users/jaina13/myPART/MethylationData/CS027417_CS027412_Reilly_210609/CS027417_Reilly_210609/"
sampleSheet<-read.table(paste0(baseDir,"CS027417_Reilly_210609.csv"),sep = ",")
#d<-list.files(paste0(baseDir,"image data/"))
#o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/SampleTPMValues.rda"))
f<-read.xlsx("~/myPART/NGS clinical data for Ashish 4-26-21v2-no_name.xlsx")
RTCancerMapping<-f[,c("Subject.Code","cancer/tumor.name.-.MyPART.Tumor.Pathologies")]
#RGset <- read.metharray.exp(paste0(baseDir,"image data/",d[1]))
#myLoad <- champ.load(paste0(baseDir,"image data/"),arraytype="EPIC")
myImport <- champ.import(paste0(baseDir,"Gabe_Data/image data/"),arraytype="EPIC")
pd<-myImport$pd
pd$RTNo<-unlist(lapply(pd$Sample_name,FUN = function(x){unlist(str_split(x,pattern = "_"))[2]}))
pdCancer<-merge(pd,RTCancerMapping,by.x="RTNo",by.y="Subject.Code",all.x=TRUE,all.y=FALSE)
pdCancer$Sample_Group<-pdCancer$`cancer/tumor.name.-.MyPART.Tumor.Pathologies`
pdCancer[is.na(pdCancer$Sample_Group),]$Sample_Group<-"Control Tumor"
myImport$pd <- pdCancer

##Filter out the control tumor samples
filtSamples<-myImport$pd$Sample_Group != "Control Tumor"
# myLoad$beta<-myLoad$beta[filtSamples,]
# myLoad$intensity<-myLoad$intensity[filtSamples,]
# myLoad$pd<-myLoad$pd[filtSamples,]
myLoad <- champ.filter(beta = myImport$beta[,filtSamples],pd = myImport$pd[filtSamples,],intensity = myImport$intensity[,filtSamples],
                       detP = myImport$detP[,filtSamples],beadcount = myImport$beadcount[,filtSamples],detPcut = 0.01,
                       SampleCutoff=0.1,arraytype="EPIC")

#champ.QC(resultsDir = paste0(baseDir,"CHAMP_QCimages"))
champ.QC(resultsDir = paste0(baseDir,"CHAMP_QCimages_group"))
myNorm <- champ.norm(beta=myLoad$beta,method="BMIQ" ,arraytype="EPIC",cores=5,resultsDir = paste0(baseDir,"CHAMP_NormImages_group"))
#champ.QC(beta = myNorm,pheno = myLoad$pd$Sample_Group,resultsDir = paste0(baseDir,"CHAMP_QCimages_norm"))
champ.QC(beta = myNorm,pheno = myLoad$pd$Sample_Group,resultsDir = paste0(baseDir,"CHAMP_QCimages_norm_group"))
champ.SVD(beta=myNorm,pd=myLoad$pd[,1:5],resultsDir = paste0(baseDir,"CHAMP_SVDImages_group/"))
#myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))

##For CNV
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group,control = FALSE,arraytype = "EPIC",resultsDir = paste0(baseDir,"CHAMP_CNVImages"))

##Mapping EPIC probes to genes
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
data("Other")
EPIC <- as.data.frame(Other)
gene<-"PTPRU"
cancer<-"Gastrointestinal stromal tumor"
geneLinkedProbesGrepl<-(EPIC[grepl(gene,EPIC$UCSC_RefGene_Name),])
g<-apply(geneLinkedProbesGrepl,1,FUN=function(x){if(gene %in% unlist(stringr::str_split(x["UCSC_RefGene_Name"],";"))){return(TRUE)}else{return(FALSE)}})
geneLinkedProbes<-intersect(row.names(geneLinkedProbesGrepl[g,]),row.names(myNorm))
#methylationData<-myNorm[geneLinkedProbes,myLoad$pd$Sample_Group == cancer]
#methylationData<-myNorm[geneLinkedProbes,]
methylationData<-myNorm[row.names(myDMP$Other_to_Gastrointestinal.stromal.tumor),]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)
n <- distinctColorPalette(length(levels(as.factor(myLoad$pd$Sample_Group))))
names(n)<-(levels(as.factor(myLoad$pd$Sample_Group)))
col_diagnosis <-n[myLoad$pd$Sample_Group]
#col_age = colorRamp2(c(0,max(features$Age,na.rm = T)), c("white","green"))
#col_rin = colorRamp2(c(0,max(features$RIN)), c("white","yellow"))
ha = HeatmapAnnotation(
  df = myLoad$pd$Sample_Group,
  col = list(Sample_Group = col_diagnosis#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             #Age = col_age,
             #Sex = c("Male"="green","Female"="brown"),
             #RIN = col_rin,
             #Race =c("White"="#B35806","Unknown"="#FDBC6B","Other"="white","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue")
             #cluster = col_cluster
  )
)
cheatmap <- ComplexHeatmap::Heatmap(methylationData,
                                    col=colorRamp2(seq(0, 1, length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=F,
                                    show_row_names = F,
                                    column_names_rot = 45,
                                    cluster_rows = T,
                                    #cluster_columns = T,
                                    column_order = order(as.factor(myLoad$pd$Sample_Group)),
                                    show_heatmap_legend = T) # Turning off to control the placement
pdf(paste0(baseDir,"MethylationHeatmap-",gene,".pdf"),width = 10,height = 7)
draw(cheatmap, show_annotation_legend = T)
dev.off()


###Differential Methylated probes analysis
DMPList<-list()
tableTypes<-table(myLoad$pd$Sample_Group)
names<-c(names(tableTypes[tableTypes != 1]))
for(name in names)
{
  diagnosis<-myLoad$pd$Sample_Group
  diagnosis[diagnosis != name]<-paste("Other")
  diagnosis[diagnosis == name]<-str_replace_all(name,pattern = " ",replacement = ".")
  #diagnosis<-as.factor(diagnosis)
  diagnosis<-as.factor(diagnosis)
  diagnosis<-relevel(diagnosis,ref = "Other")
  myDMP <- champ.DMP(beta = myNorm,pheno=diagnosis,arraytype = "EPIC")
  DMPList[[name]]<-myDMP
}
save(DMPList,myLoad,myNorm,file = paste0(baseDir,"DMPList.rda"))

###Map probes to genes
fcThres<-log2(1.5)
d<-DMPList$Neuroendocrine$Other_to_Neuroendocrine
g<-EPIC[row.names(d[d$logFC >= fcThres,]),]
genesHyperMethylated<-unique(unlist(apply(g,1,FUN=function(x){return(unlist(stringr::str_split(x["UCSC_RefGene_Name"],";")))})))
g<-getGOEnrichmentWebGestaltRWGCNA(genesHyperMethylated,"geneontology_Biological_Process_noRedundant","NET-HyperMeth",baseDir)
g<-getGOEnrichmentWebGestaltRWGCNA(genesHyperMethylated,"pathway_Reactome","NET-HyperMeth",baseDir)
g<-EPIC[row.names(d[d$logFC <= -fcThres,]),]
genesHypoMethylated<-unique(unlist(apply(g,1,FUN=function(x){return(unlist(stringr::str_split(x["UCSC_RefGene_Name"],";")))})))
g<-getGOEnrichmentWebGestaltRWGCNA(genes,"geneontology_Biological_Process_noRedundant","NET-HypoMeth",baseDir)
g<-getGOEnrichmentWebGestaltRWGCNA(genes,"pathway_Reactome","NET-HypoMeth",baseDir)


myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=NULL, arraytype="EPIC",adjPval=0.05, method="gometh")

