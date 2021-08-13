library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)

###Overlap between Starfusion and Arriba results
g1<-getSheetNames("/Users/jaina13/myPART/FusionAnalysis/CombinedFusionGenesWONormal/unique.fusion.genes.WONormal.xlsx")
g2<-getSheetNames("/Users/jaina13/myPART/FusionAnalysis/arriba_results/unique.fusion.protein-coding.genes.arriba.xlsx")
l2<-list()
for(name in intersect(g1,g2))
{
  #name<-str_replace_all(name,pattern = "[ /]",replacement = ".")
  starfusion<-read.xlsx("/Users/jaina13/myPART/FusionAnalysis/CombinedFusionGenesWONormal/unique.fusion.genes.WONormal.xlsx",sheet = name)
  starfusion$Fusion<-paste0(starfusion$Gene1,"-",starfusion$Gene2,"-",starfusion$Gene1_Chr,":",starfusion$Gene1_Position,"-",starfusion$Gene2_Chr,":",starfusion$Gene2_Position)
  arriba<-read.xlsx("/Users/jaina13/myPART/FusionAnalysis/arriba_results/unique.fusion.protein-coding.genes.arriba.xlsx",sheet = name)
  m<-merge(starfusion,arriba,by="Fusion")
  if(nrow(m)>0)
  {
    l2[[name]]<-m
    mergeFusionDF<-m
    curr_fusions<-mergeFusionDF[!duplicated(mergeFusionDF$Fusion),]
    left_break <- cbind(curr_fusions$Gene1_Chr,curr_fusions$Gene1_Position)
    right_break <- cbind(curr_fusions$Gene2_Chr,curr_fusions$Gene2_Position)
    gene_names <- cbind(curr_fusions$Gene1,curr_fusions$Gene2)
    sampleid <- (curr_fusions$Sample)
    #SpanningFragCount <- curr_fusions$SpanningFragCount.x
    #JunctionReadCount <- curr_fusions$JunctionReadCount.x
    #ffpm <- curr_fusions$FFPM.x
    nos <- curr_fusions$Number
    TCGAAnnot<-curr_fusions$TCGAAnnot
    CCLEAnnot<-curr_fusions$CCLEAnnot
    if(nrow(curr_fusions) >1)
    {
      return_df <- as.matrix(cbind(gene_names,left_break[,1:2], right_break[,1:2],nos,sampleid,TCGAAnnot,CCLEAnnot))
      colnames(return_df) <- c("Gene1","Gene2","Gene1_Chr","Gene1_Position","Gene2_Chr","Gene2_Position","No. Of Samples","sample_id","TCGAAnnot","CCLEAnnot")
    }else
    {
      return_df <- t(as.matrix(c(gene_names,left_break[,1:2], right_break[,1:2],nos,sampleid,TCGAAnnot,CCLEAnnot)))
      colnames(return_df) <- c("Gene1","Gene2","Gene1_Chr","Gene1_Position","Gene2_Chr","Gene2_Position","No. Of Samples","sample_id","TCGAAnnot","CCLEAnnot")
    }
    links_data<-(data.frame(return_df))
    lwd<-as.numeric(links_data$No..Of.Samples)*0.8
    col<-apply(links_data,1,FUN=function(x){
      if(x["TCGAAnnot"] == "Y" && x["CCLEAnnot"] == "Y"){return("pink")}
      if(x["TCGAAnnot"] == "Y") {return("red")} else if(x["CCLEAnnot"] == "Y") {return("blue")} else {return("black")}
    })

    bed1 <- data.frame(chr=links_data$Gene1_Chr,
                       start=as.numeric(links_data$Gene1_Position),
                       end=as.numeric(links_data$Gene1_Position)+1

    )
    #write.table(links_data,file = paste0(mywd,"UniqueFusions-",name,"-Gene1.bed"),sep = "\t",quote = F,row.names = F)
    bed2 <- data.frame(chr=links_data$Gene2_Chr,
                       start=as.numeric(links_data$Gene2_Position),
                       end=as.numeric(links_data$Gene2_Position)+1
    )
    #write.table(links_data,file = paste0(mywd,"UniqueFusions-",name,"-Gene2.bed"),sep = "\t",quote = F,row.names = F)

    names_bed <- rbind(unique(cbind(links_data$Gene1_Chr, links_data$Gene1_Position, as.numeric(links_data$Gene1_Position)+1, links_data$Gene1)),
                       unique(cbind(links_data$Gene2_Chr, links_data$Gene2_Position, as.numeric(links_data$Gene2_Position)+1, links_data$Gene2)))
    names_bed <- data.frame(names_bed, stringsAsFactors = F)
    names_bed[,2] <- as.numeric(names_bed[,2])
    names_bed[,3] <- as.numeric(names_bed[,3])
    names_bed <- names_bed[!duplicated(names_bed[,4]),]

    #cluster_colors = c(MO1="#4DAF4A",MO2="#377EB8",MO3="orange",MO4="#E41A1C","NA"="grey70")
    mychrs=sort(unique(c(links_data$Gene1_Chr,links_data$Gene2_Chr)))
    mychrs=mychrs[order(nchar(mychrs), mychrs)]

    fileName<-paste0("../fusions.common.star.arriba",name,".pcfiltered.pdf")
    pdf(paste0(mywd,fileName),width=8, height=8)
    circos.par(start.degree = 91,gap.degree=c(rep(1, length(mychrs)-1), 8))
    circos.initializeWithIdeogram(chromosome.index = mychrs,species = "hg38")
    circos.genomicLabels(names_bed, labels.column = 4, side = "outside",
                         col = "grey30", line_col = "grey50")
    circos.genomicLink(bed1, bed2, lwd=lwd,col = col)
    circos.clear()
    title(paste0("Fusions in Starfusion and Arriba"))
    dev.off()
  }
}

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(l2, paste0("/Users/jaina13/myPART/FusionAnalysis/unique.fusion.protein-coding.genes.arriba-starfusion.combined.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)


"create_report test.maflite.tsv hg19.fasta --flanking 20 --sequence 1 --begin 2 --end 3 --info-columns chr start end ref_allele alt_allele --tracks
Tumor_Match_Run_Pipeliner_1/Sample_10_HuDTumMix1_1_5_S10.recal.bam Tumor_Match_Run_Pipeliner_1/Sample_44_HuDNA6964PB1_1_S44.recal.bam
Tumor_Match_Run_Pipeliner_1/Sample_30_HuDTumMix1_1_1_S30.recal.bam
Tumor_Match_Run_Pipeliner_1/Sample_64_HuDNA6964PB1_2_S64.recal.bam
--output example_tab.html"

"create_report Control.SNVs.txt /data/tandonm/pl_test_data/human/genome/hg38/Homo_sapiens_assembly38.fasta --flanking 25 --sequence 1 --begin 2 --end 3 --info-columns chromosome start end TranscriptChange Gene CellLine --tracks Tumor_Match_Run_Pipeliner_1/Sample_10_HuDTumMix1_1_5_S10.recal.bam Tumor_Match_Run_Pipeliner_1/Sample_30_HuDTumMix1_1_1_S30.recal.bam Tumor_Match_Run_Pipeliner_1/Sample_44_HuDNA6964PB1_1_S44.recal.bam Tumor_Match_Run_Pipeliner_1/Sample_64_HuDNA6964PB1_2_S64.recal.bam --output control_wes_igv_report_tab.html"
