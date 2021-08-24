library(WGCNA)

################Supplementary Figure################
####################Module Heatmap######################################
library(WGCNA)
# The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
lnames = load(file = paste0(parent,"MyPART-dataInput.RData"));
lnames = load(file = paste0(parent,"MyPART-networkConstruction-merged.RData"));
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

##Renaming the module Names
moduleNames<-colnames(MEs)
names(moduleNames)<-paste0("M",c(1:ncol(MEs)))
colnames(MEs)<-names(moduleNames[moduleNames == colnames(MEs)])
colorNames<-paste0("M",c(1:ncol(MEs)))
names(colorNames)<-moduleNames

datTraits<-datTraits[row.names(datExpr),]
moduleTraitCor = cor(MEs, datTraits[,c(5,6,13:14,10:12,15:16)], use = "p");#c(5:6,10:12)
#moduleTraitCor = cor(MEs, datTraits[,c(5:8,14)], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# adjmoduleTraitPvalue = apply(moduleTraitPvalue,2,function(x){
#   p.adjust(x,method = "BH")
# })
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#png(filename=paste0(parent,"moduleSignificanceHeatmap.png"), width = 4500 , height = 3500,units = "px",res = 300)
pdf(paste0(parent,"moduleSignificanceHeatmap.pdf"),height = 10,width = 9)
WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                      xLabels = names(datTraits[,c(5,6,13:14,10:12,15:16)]),
                      #yLabels = modulesName,
                      yLabels = names(MEs),
                      ySymbols = names(MEs),
                      colorLabels = FALSE,
                      #colors = colorRampPalette(c("white","red"),space = "Lab")(50),
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix,
                      setStdMargins = TRUE,
                      cex.text = 0.7,
                      cex.lab = 1,
                      xLabelsAngle = 45,
                      zlim = c(-1,1),
                      main = paste("Module-Trait Significance"))
dev.off()



###########Figure 2B##################
parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/"
l<-load(file = paste0(parent,"MyPART-networkConstruction-merged.RData"))
l1<-load(file = paste0(parent,"MyPART-dataInput.RData"))
datTraits<-datTraits[row.names(datExpr),]
datTraits<-datTraits[order(datTraits$Diagnosis),]
datExpr<-datExpr[row.names(datTraits),]
MEs<-MEs[row.names(datTraits),]
colnames(MEs)<-names(moduleNames[moduleNames == colnames(MEs)])
colorh1 <-colorNames[paste0("ME",moduleColors)]

##Module heatmap and eigengene
library(randomcoloR)
colors <- c("#80D9D8","#D98B69","#B7DE9C","#D6CDCD","#AFE456","#C457D3","#B096D3")#distinctColorPalette(length(unique(datTraits$Diagnosis)))#c("orange","tan")#
names(colors)<-unique(datTraits$Diagnosis)
###Generate a legend for type of samples in the data###########
png(filename=paste0(parent,"eigenGeneLegend-Diagnosis.png"), width = 1600 , height = 1600,units = "px",res = 300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(colors), pch=15, pt.cex=3, cex=1.5, bty='n',col = colors)
mtext("Diagnosis", at=0.2, cex=2)
dev.off()
#################################################################
groupColors<-colors[as.character(datTraits$Diagnosis)]
#colo<-c("brown","cyan","red","turquoise","yellow","salmon")#unique(colorh1)##c("grey60")#
colo<-c("M1","M3","M4","M5","M15","M16")
#datTraits$RTNo<-unlist(lapply(row.names(datTraits), function(x){unlist(str_split(x,pattern = "_"))[2]}))
#datTraits$RTNo[datTraits$Diagnosis == "Normal"]<-datTraits$Tissue[datTraits$Diagnosis == "Normal"]
for (c in colo)
{
  which.module=c
  #moduleNumber<-"M14"
  #ME=MEs[row.names(datTraits), paste0("ME",which.module)]
  ME=MEs[row.names(datTraits), which.module]
  #names(ME)<-datTraits$RTNo#unlist(lapply(row.names(datTraits), function(x){unlist(str_split(x,pattern = "_"))[2]}))
  #names(groupColors)<-names(ME)
  #ME<-sort(ME,decreasing = T)
  #groupColors<-groupColors[names(ME)]
  #png(filename=paste(parent,"M",which.module,"-ModulesExpression.png",sep = ""), width = 2500 , height = 1600,units = "px",res = 300)
  pdf(file =paste(parent,"M",which.module,"-ModulesExpression.pdf",sep = ""), width = 8.39 , height = 5.33)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  #WGCNA::plotMat(t(scale(datExpr[,colorh1==which.module]) ),nrgcols=30,rlabels=F,rcols=which.module,main=paste0("Module ",which.module), cex.main=2)
  plotMat(t(scale(datExpr[,colorh1==which.module]) ),nrgcols=30,rlabels=F,rcols=names(colorNames[which.module]),main=paste0("Module ",which.module), cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=groupColors, main="", cex.main=2,ylab="Eigengene Expression",xlab="",las = 2,cex.names = 0.7)
  #barplot(ME, col="grey", main="", cex.main=2,ylab="Eigengene Expression",xlab="",las = 2,cex.names = 0.7)
  dev.off()
}

plotMat<-function(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="", ...)
{
  n<-nrow(x)
  p<-ncol(x)
  image(1:p,1:n,t(x[n:1,]),col=colorRampPalette(c("#E9C226","#000099"))(nrgcols),axes=FALSE, xlab="", ylab="", ... )

  if(length(ccols)==1){
    axis(3,at=1:p,labels=clabels,las=2,cex.axis=0.6,col.axis=ccols)
  }

  if(length(ccols)==p){
    cols<-unique(ccols)
    for(i in 1:length(cols)){
      which<-(1:p)[ccols==cols[i]]
      axis(3,at=which,labels=clabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
    }
  }

  if(length(rcols)==1){
    axis(2,at=n:1,labels=rlabels,las=2,cex.axis=0.6,col.axis=rcols)
  }

  if(length(rcols)==n){
    cols<-unique(rcols)
    for(i in 1:length(cols)){
      which<-(1:n)[rcols==cols[i]]
      axis(2,at=(n:1)[which],labels=rlabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
    }
  }
  mtext(title,side=3,line=3)
  box()
}


###########Figure 2A##################
###########GO Ontology Terms###########
modulesGOTerms<-c("brown","yellow","salmon","red")
names(modulesGOTerms)<-c("M4","M5","M15","M16")
#####Read top 5 GO terms of each of the four modules
top5TermsAll<-list()
for(co in modulesGOTerms)
{
  data<-read.table(paste0(parent,"EnrichmentAnalysis/Project_",co,"_geneontology_Biological_Process_noRedundant/enrichment_results_",co,"_geneontology_Biological_Process_noRedundant.txt"),sep = "\t",header = T)
  enrichmentRes<-data[,c(1,2,7,9)] %>% dplyr::filter(FDR <= 0.05)
  top5terms<-enrichmentRes[order(enrichmentRes$enrichmentRatio,decreasing = TRUE),][1:5,]
  top5terms$Module<-names(modulesGOTerms[modulesGOTerms == co])
  top5TermsAll[[co]]<-top5terms
}
top5TermsAllDF<-do.call("rbind",top5TermsAll)
top5TermsAllDF$description[6]<-paste0(top5TermsAllDF$description[6]," ")
top5TermsAllDF$description <- tools::toTitleCase(top5TermsAllDF$description)
#data[6,"V3"]<-data[6,"V3"]-5
#data[7,"V3"]<-data[7,"V3"]-2
p<-ggplot(top5TermsAllDF,aes(x=reorder(description,-c(1:20)),y=enrichmentRatio,fill=Module))+
  geom_bar(stat = 'identity', position = position_dodge(width=1))+
  coord_flip()+
  theme_bw()+
  # geom_hline(yintercept = 3.9,linetype="solid") +
  # geom_hline(yintercept = 4.1,linetype="solid") +
  # geom_hline(yintercept = 5.9,linetype="solid") +
  # geom_hline(yintercept = 6.1,linetype="solid") +
  #guides(fill = "none")+
  #scale_y_continuous(breaks = 1:8, labels = c(1:3,"Break",7,"Break",12:13)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme(text=element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20),axis.text = element_text(size=18,colour = "black"),axis.title = element_text(size=20,colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #labs(x='', y = bquote(-Log[10]~'(FDR)'))
  labs(x='', y = "Enrichment Ratio")
  #facet_wrap(~Module,ncol = 4)
#ggtitle('Module Top 5 GO Biological Process terms')
ggsave(paste0(parent,"Modules.GOterms.pdf"),width = 6.6,height = 3.7,units = "in",scale = 2,p)


################Supplementary Figures################
################Soft Threshold Plot#################
lnames <- load(file = paste0(parent,"MyPART-dataInput.RData"))
typeOfCorrelation<-"signed hybrid"
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = typeOfCorrelation,corFnc = "bicor",corOptions = list(maxPOutliers =0.1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = typeOfCorrelation)
# Plot the results:
pdf(paste0(parent,"softThresholdplot.pdf"),height = 4.7,width = 8.2)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.7;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

####################Module Distance Dendrogram######################################
# Calculate dissimilarity of module eigengenes and keeps track of sign of correlation.
dissimME=(1-t(cor(MEs, method="p")))/2
#dissimME[dissimME < 0] <- 0
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
pdf(paste0(parent,"moduleClusteringDendogram.pdf"),height = 5.5,width = 8.2)
plot(hclustdatME, main="Clustering of modules based on eigengenes")
dev.off()


#######################Pathway Enrichment Dotplots###################################
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

data_folder = "/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/EnrichmentAnalysis/"
gene_sig_cutoff = 0.05

topPathways<-10
#modules<-c("lightcyan","turquoise","blue","green")
modulesPathTerms<-c("brown","yellow","red","turquoise","cyan")
names(modulesPathTerms)<-c("M4","M5","M16","M3","M1")
#modules<-c("turquoise","cyan","brown","yellow","red")
all_res <- lapply(modulesPathTerms, FUN=function(x){
  t<-readr::read_tsv(file = paste0(data_folder,"/Project_",x,"_pathway_Reactome/enrichment_results_",x,"_pathway_Reactome.txt"),progress = FALSE)
  t$contrast<-names(modulesPathTerms[modulesPathTerms == x])
  t<-t[order(t$FDR),]
  if(nrow(t)>=topPathways)
  {
    t<-t[1:topPathways,]
  }
  return(t)
})

names(all_res)<-names(modulesPathTerms)
#names(all_res) <- lapply(res_files, function(x) {unlist(lapply(strsplit(basename(x),"\\."),"[[",1))})

#con<-factor(enrichment_input$contrast,levels=c("M4","M5","M16"))
#enrichment_input$contrast<-con
source("/Users/jaina13/myPART/MyPART-Analysis/RScripts/helper_functions.enrichment.R")
enrichment_input <- do.call(rbind,all_res)
plot_all.mod <- add_category_to_dotplot_Webgestalt(enrichment_input, show_n_path = 100)
#plot_all.mod <- add_category_to_dotplot(ck.reactome.all, show_n_path = 100)
# plot_all.mod <- add_category_to_dotplot(ck.reactome.all, show_n_path = 20)

pdf(file = paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/pathways_by_cluster.both.category.top10.pdf"),width = 24, height=15)
# pdf(file = paste0("pathways_by_cluster.both.category.top_20.pdf"),width = 16, height=6)
myplot_with_cat <- plot_all.mod[[2]]
myplot_with_cat <- myplot_with_cat +
  #scale_fill_manual(values=c(up="green3",down="red")) +
  theme_bw() +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme(text=element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20),axis.text = element_text(size=15,colour = "black"),axis.title = element_text(size=20,colour = "black"))
  #theme(axis.text.x = element_text(angle = 20, hjust = 1))
print(myplot_with_cat)
dev.off()

#######################GO Enrichment TreeMap###################################
modulesGOTerms<-c("brown","yellow","salmon","red","turquoise","cyan")
names(modulesGOTerms)<-c("M4","M5","M15","M16","M3","M1")
topPathways<-20
all_res <- lapply(modulesGOTerms, FUN=function(x){
  t<-readr::read_tsv(file = paste0(data_folder,"/Project_",x,"_geneontology_Biological_Process_noRedundant/enrichment_results_",x,"_geneontology_Biological_Process_noRedundant.txt"),progress = FALSE)
  t$contrast<-names(modulesPathTerms[modulesGOTerms == x])
  t<-t[order(t$FDR),]
  if(nrow(t)>=topPathways)
  {
    t<-t[1:topPathways,]
  }
  return(t)
})
names(all_res)<-names(modulesGOTerms)

library(rrvgo)
go_analysis <- all_res$M1 #read.delim(system.file("extdata/example.txt", package="rrvgo"))
#go_analysis$description<-tools::toTitleCase(go_analysis$description)
simMatrix <- calculateSimMatrix(go_analysis$geneSet,orgdb="org.Hs.eg.db",ont="BP",method="Rel")
minValue<-ifelse(min(go_analysis$FDR[go_analysis$FDR !=0]) < 10e-7,min(go_analysis$FDR[go_analysis$FDR !=0]),10e-7)
scores <- setNames(-log10(go_analysis$FDR+(minValue)), go_analysis$geneSet)
reducedTerms <- reduceSimMatrix(simMatrix,scores,threshold=0.7,orgdb="org.Hs.eg.db")
reducedTerms$term<-tools::toTitleCase(reducedTerms$term)
reducedTerms$parentTerm<-tools::toTitleCase(reducedTerms$parentTerm)
pdf(file = paste0(parent,"GOTerms.Treemap.M1.pdf"),width = 10, height=10)
treemapPlot(reducedTerms,fontsize.labels = 20)
dev.off()

