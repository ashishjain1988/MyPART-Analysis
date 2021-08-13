library(WGCNA)
source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helperFunctions-WGCNA.R")

typeOfCorrelation<-"signed hybrid"
parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/"
l<-load(file = paste0(parent,"../vsdNormalizedCounts.rda"))
datExpr0 <- vsdNormalizedCounts

rnaSamples<-row.names(features[features$Tissue == "Human Universal Reference Total RNA",])
###Removing the samples with bad quality (low TIN) from the analysis
lowQualitySamples<-intersect(row.names(features),c("47_RT00085_X007R_resent","76_RT00105_X336R"))
datExpr0 <- datExpr0[,!(colnames(datExpr0) %in% c(rnaSamples,lowQualitySamples))]

# filtSamples<-datTraits$Diagnosis %in% c("Adrenocortical carcinoma","Gastrointestinal stromal tumor","Carcinoid tumor","Normal")
# datExpr0<-datExpr0[,filtSamples]
# datTraits<-datTraits[colnames(datExpr0),]
# datExpr0<-datExpr0[apply(datExpr0,MARGIN = 1, function(x) any(x >= 1)), ]

#=====================================================================================
#
#  Code chunk: Processing the gene expression data before analysis.
#
#=====================================================================================
variance<-apply(datExpr0,MARGIN = 1, function(x) var(x))
datExpr0<-datExpr0[names(variance[variance>quantile(variance,0.75)]), ]

##For my Data
datExpr<-as.data.frame(t(as.matrix(datExpr0)))
##Loading Clinical Data
datTraits <- read.csv(paste0(parent,"traits.txt"),row.names = 1,header = T,sep = "\t");
datTraits<-datTraits[row.names(datExpr),]
save(datExpr,datTraits, file = paste0(parent,"MyPART-dataInput.RData"))
#=====================================================================================
#
#  Code chunk: Calculating the scale free thrshold for the co-expression network.
#
#=====================================================================================
##Create the co-expression network
# Load the WGCNA package
library(WGCNA)
lnames <- load(file = paste0(parent,"MyPART-dataInput.RData"))
powers = c(c(1:10), seq(from = 12, to=25, by=2))
#registerDoSEQ()
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = typeOfCorrelation)
png(filename=paste(parent,"scalefreethreshold.png",sep = "/"), width = 3500 , height = 2000,units = "px",res = 300)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.9,col="red")
abline(h=0.8,col="red")
#abline(h=0.7,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

#=====================================================================================
#
#  Code chunk: Construction of co-expression network.
#
#=====================================================================================
##Check for the power from the power and mean connectivity plot
softPower = 5
adjacency = WGCNA::adjacency(datExpr, type = typeOfCorrelation, power = softPower);
TOM = TOMsimilarity(adjacency,TOMType = "signed");
colnames(TOM) = rownames(TOM) =rownames(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
minModuleSize = 50;
# Module identification using dynamic tree cut using dissimilairt based on TOM
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
module_colors= names(table(dynamicColors))
##Writing the unmerged modules in file having prfic module followed by module color.
# SubGeneNames<-colnames(datExpr)
# module_colors= names(table(dynamicColors))
# for (color in module_colors){
#   module=SubGeneNames[which(dynamicColors==color)]
#   write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
# }
# # Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
save(MEs, dynamicColors, geneTree, file = paste0(parent,"/MyPART-networkConstruction-unmerged.RData"))
#=====================================================================================
#
#  Code chunk: Merging of Module based on correlation between them.
#
#=====================================================================================
l1<-load(file = paste0(parent,"MyPART-dataInput.RData"))
load(file = paste0(parent,"/MyPART-networkConstruction-unmerged.RData"))
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#sizeGrWindow(12, 9)
pdf(file = paste0(parent,"geneDendro-3.pdf"), wi = 9, he = 6)
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(geneTree, mergedColors,
                    "Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
##Writing the Merged modules in file having prfic module followed by module color.
SubGeneNames<-colnames(datExpr)
module_colors= names(table(mergedColors))
for (color in module_colors){
  module=SubGeneNames[which(mergedColors==color)]
  print(color)
  n<-extractGeneNames(module)
  write.table(n$GeneName, paste(parent,"module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
  g<-getGOEnrichmentWebGestaltRWGCNA(n$GeneName,"geneontology_Biological_Process_noRedundant",color,parent)
  g<- getGOEnrichmentWebGestaltRWGCNA(n$GeneName,"pathway_Reactome",color,parent)
}
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
save(MEs, moduleColors, geneTree, file = paste0(parent,"MyPART-networkConstruction-merged.RData"))

#=====================================================================================
#
#  Code chunk: Module analysis based on eigengene.
#
#=====================================================================================
l<-load(file = paste0(parent,"MyPART-networkConstruction-merged.RData"))
l1<-load(file = paste0(parent,"MyPART-dataInput.RData"))
datTraits<-datTraits[row.names(datExpr),]
datTraits<-datTraits[order(datTraits$Diagnosis),]
datExpr<-datExpr[row.names(datTraits),]
MEs<-MEs[row.names(datTraits),]
colorh1 <-moduleColors
#MEsVariance<- MEList$varExplained
#names(MEsVariance)<-colnames(MEs)
# Calculate dissimilarity of module eigengenes and keeps track of sign of correlation.
dissimME=(1-t(cor(MEs, method="p")))/2
#dissimME[dissimME < 0] <- 0

# Plot the eigengene dendrogram
hclustdatME=hclust(as.dist(dissimME), method="average" )
png(filename=paste(parent,"moduleDendogram.png",sep = "/"), width = 3000 , height = 2000,units = "px",res = 300)
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering of modules based on eigengenes",cex = 1.5,cex.main=1.5,cex.lab=1.4)
dev.off()

##Pairwise Scatter plot of the samples along the module eigengenes
# sizeGrWindow(8,9)
# plotMEpairs(MEs)
#plotMEpairs(datME,datTraits$ESCd)

##Module heatmap and eigengene
#sizeGrWindow(8,7);
library(randomcoloR)
colors <- distinctColorPalette(length(unique(datTraits$Diagnosis)))#c("orange","tan")#
names(colors)<-unique(datTraits$Diagnosis)
groupColors<-colors[as.character(datTraits$Diagnosis)]
colo<-unique(colorh1)#c("black","blue","brown","pink","yellow","purple","royalblue")#c("grey60")#
datTraits$RTNo<-unlist(lapply(row.names(datTraits), function(x){unlist(str_split(x,pattern = "_"))[2]}))
datTraits$RTNo[datTraits$Diagnosis == "Normal"]<-datTraits$Tissue[datTraits$Diagnosis == "Normal"]
for (c in colo)
{
  which.module=c
  #moduleNumber<-"M14"
  ME=MEs[row.names(datTraits), paste0("ME",which.module)]
  names(ME)<-datTraits$RTNo#unlist(lapply(row.names(datTraits), function(x){unlist(str_split(x,pattern = "_"))[2]}))
  names(groupColors)<-names(ME)
  #ME<-sort(ME,decreasing = T)
  #groupColors<-groupColors[names(ME)]
  png(filename=paste(parent,"M",which.module,"-ModulesExpression.png",sep = ""), width = 2500 , height = 1600,units = "px",res = 300)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[,colorh1==which.module]) ),nrgcols=30,rlabels=F,rcols=which.module,main=paste0("Module ",which.module), cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=groupColors, main="", cex.main=2,ylab="Eigengene Expression",xlab="",las = 2,cex.names = 0.7)
  #barplot(ME, col="grey", main="", cex.main=2,ylab="Eigengene Expression",xlab="",las = 2,cex.names = 0.7)
  # abline(v=19.2)
  # abline(v=(19.2+8.6))
  # abline(v=(19.2+8.7+16.5))
  dev.off()
}
#=====================================================================================
#
#  Code chunk: Assigning genes to multiple modules based on the correlation with
#              module eigengene.
#
#=====================================================================================
##Group genes based on KME values
datKME<-signedKME(datExpr, MEs, outputColumnName="MM.")
for (color in colnames(datKME)){
  MMFilt<-(datKME[(datKME[,color]) >= 0.75,which(colnames(datKME)==color),drop=FALSE])
  n<-extractGeneNames(row.names(MMFilt))
  write.table(n$GeneName, paste(parent,"datKME/","module_",color, ".txt",sep=""), sep="\t", row.names = FALSE, col.names=FALSE,quote=FALSE)
  #write.table(n$GeneName, paste(parent,"module_",color, ".GREAT.txt",sep=""), sep="\t", row.names = FALSE, col.names=FALSE,quote=FALSE)
  g<-getGOEnrichmentWebGestaltRWGCNA(n$GeneName,"geneontology_Biological_Process_noRedundant",color,paste0(parent,"datKME/"))
  g<- getGOEnrichmentWebGestaltRWGCNA(n$GeneName,"pathway_Reactome",color,paste0(parent,"datKME/"))
}

#=====================================================================================
#
#  Code chunk: Linking module to cell type based on correlation between cell-type data
#              gene expression and module eigengene.
#
#=====================================================================================
#parent = "/Users/jaina13/myPART/count-file/WGCNA/";
#setwd(workingDir);
library(WGCNA)
# The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
lnames = load(file = paste0(parent,"MyPART-dataInput.RData"));
lnames = load(file = paste0(parent,"MyPART-networkConstruction-merged.RData"));
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
##For unmerged dataset
#moduleColors<-dynamicColors
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#datTraits = read.csv(paste0(parent,"traits.txt"),row.names = 1,header = T,sep = "\t")
datTraits<-datTraits[row.names(datExpr),]
#datTraits$AdrenalGland<-as.numeric(datTraits$Tissue == "adrenal")
moduleTraitCor = cor(MEs, datTraits[,c(5,6,13,10:12,14)], use = "p");#c(5:6,10:12)
#moduleTraitCor = cor(MEs, datTraits[,c(5:8,14)], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# adjmoduleTraitPvalue = apply(moduleTraitPvalue,2,function(x){
#   p.adjust(x,method = "BH")
# })
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
png(filename=paste0(parent,"moduleSignificanceHeatmap.png"), width = 4500 , height = 3500,units = "px",res = 300)
WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                      xLabels = names(datTraits[,c(5,6,13,10:12,14)]),
                      #yLabels = modulesName,
                      yLabels = names(MEs),
                      ySymbols = names(MEs),
                      colorLabels = FALSE,
                      #colors = colorRampPalette(c("white","red"),space = "Lab")(50),
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix,
                      setStdMargins = TRUE,
                      cex.text = 0.9,
                      cex.lab = 1,
                      xLabelsAngle = 45,
                      zlim = c(-1,1),
                      main = paste("Module-Traits Significance"))
dev.off()

#=====================================================================================
#
#  Code chunk: Clustering samples based on the eigengenes across the modules
#
#=====================================================================================
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)

hm_data<-as.matrix(t(MEs[row.names(featuresFilt),]))
column_annotations = HeatmapAnnotation(df = featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","green"))
#col_rin = colorRamp2(c(0,max(featuresFilt$RIN)), c("white","yellow"))
n <- distinctColorPalette(length(levels(as.factor(featuresFilt$Tissue))))
names(n)<-(levels(as.factor(featuresFilt$Tissue)))
col_tissue <-n[featuresFilt$Tissue]

n1 <- distinctColorPalette(length(levels(as.factor(featuresFilt$New.Diagnosis))))
names(n1)<-(levels(as.factor(featuresFilt$New.Diagnosis)))
col_disease <-n1[featuresFilt$New.Diagnosis]

df<-featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")]
colnames(df)<-c("Tissue","Diagnosis","Age","Sex","Race","TumorSite")
#col_tissue <- distinctColorPalette()
ha = HeatmapAnnotation(
  #df = features[,c(5,6,7,13,15)],
  #df = featuresFilt[,c(3,4,5,6,13,7)],
  df = df,
  col = list(Tissue = col_tissue,
             Diagnosis = col_disease,#c("Normal" = "black","Adrenocortical carcinoma" = "red"),#col_cluster,#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             Age = col_age,
             Sex = c("Male"="green","Female"="brown"),
             #RIN = col_rin,
             Race =c("White"="#B35806","Unknown"="#FDBC6B","Other"="white","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue"),
             TumorSite=c("metastasis"="green","primary site"="brown","recurrence"="white")
             #cluster = col_cluster
  )
)

cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(hm_data), na.rm = T), max(abs(hm_data), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=F,
                                    column_names_rot = 45,
                                    cluster_rows = FALSE,
                                    show_heatmap_legend = F) # Turning off to control the placement
pdf(paste0(parent,"PCsSampleClustering-Subtypes-All-L.pdf"))
draw(cheatmap, show_annotation_legend = TRUE)
dev.off()

#=====================================================================================
#
#  Code chunk: Combining genes in modules with module membership and trait
#  significance score.
#
#=====================================================================================


# l<-load(file = paste0(parent,"MyPART-networkConstruction-merged.RData"))
# l1<-load(file = paste0(parent,"MyPART-dataInput.RData"))
# datTraits<-datTraits[row.names(datExpr),]
# CancerType = as.data.frame(datTraits$NET);
# names(CancerType) = "NET"
# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# # Define numbers of genes and samples
# nGenes = ncol(datExpr);
# nSamples = nrow(datExpr);
#
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
#
# geneTraitSignificance = as.data.frame(cor(datExpr, CancerType, use = "p"))
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# names(geneTraitSignificance) = paste("GS.", names(CancerType), sep="")
# names(GSPvalue) = paste("p.GS.", names(CancerType), sep="")
#
# ##Plotting them
# module = "turquoise"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# geneNames<-colnames(datExpr)[moduleGenes]
# png(filename=paste0(parent,module,"ModuleScatterPlot.png"), width = 4500 , height = 3500,units = "px",res = 300)
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for ACC",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# dev.off()
#
# ###Values for genes in modules
# d<-data.frame(MM=geneModuleMembership[moduleGenes, column],pMM=MMPvalue[moduleGenes, column],
#               CancerType=geneTraitSignificance[moduleGenes, 1],pCancerType=GSPvalue[moduleGenes, 1],row.names = geneNames)
# d1<-d %>% filter(pMM <= 0.05 & pCancerType <= 0.05 & abs(MM) >=0.75 & abs(CancerType) >=0.5)
# n<-extractGeneNames(row.names(d1))
# d1$GeneName<-n$GeneName
# d1$EnsemblID<-n$EnsemblID
#
#
# teGenes<-getTissueSpecificGenes(rdaFilePath = "~/myPART/TissueEnrichCombineExpression.rda")
# hpa<-teGenes$HPA[,c(1:3)]
# hpa1<-merge(x=d1,y=hpa,by.x="EnsemblID",by.y="Gene",all.x=TRUE,all.y=FALSE)
# colnames(hpa1)<-c(colnames(hpa1)[1:6],"HPA-Tissue","HPA-TE-Type")
# gtex<-teGenes$GTEX[,c(1:3)]
# hpa1<-merge(x=hpa1,y=gtex,by.x="EnsemblID",by.y="Gene",all.x=TRUE,all.y=FALSE)
# colnames(hpa1)<-c(colnames(hpa1)[1:8],"GTEx-Tissue","GTEx-TE-Type")
# row.names(hpa1)<-row.names(d1)
# hpa1<-hpa1[,c(2,3,4,5,7,9)]
