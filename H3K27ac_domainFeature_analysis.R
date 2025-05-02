#!/usr/bin/env Rscript

# 3DSIM H3K27ac Domain Feature Analysis
# Author: Shengen Shawn Hu
# Date: May 2025
#
# Description:
#   This script performs downstream statistical analysis and visualization of H3K27ac-domain-level
#   features extracted from 3D-SIM images. It compares treatment effects across MCF7 cell lines
#   (WT and YS mutant) under various hormone and antagonist conditions (E2, ED, TAM, FULV).
#

library(gplots)
library(ggplot2)
library(tidyr)
library(gridExtra)

# -------------------------
# Read in collected features
# -------------------------
# features were summarized using 3DSIM_pipeline.py for each samples
domain_use<-read.table("H3K27ac_domain_features.txt")
# fetch meta data
meta_mat <- cbind(do.call(rbind, strsplit(domain_use[,"name"], "_")),domain_use[,"cellID"])
meta_mat_use <- cbind(paste(meta_mat[,1],meta_mat[,2],meta_mat[,3],meta_mat[,4],sep="_"), meta_mat)

# ---------------------------------
# Feature Analysis (batch-combined)
# ---------------------------------
# material for figure F1C, SF2
useFeatures <- c("spher","corH3K27acER","volumn","distBound")
allPV <- c() 
for(thisFeature in useFeatures ){
  usedata<- as.data.frame( cbind(meta_mat_use, domain_use[,thisFeature]) )
  colnames(usedata) <- c("name","celltype","chan1","chan2","batch","treatment","rep","cell","feature")
  usedata$feature <- as.numeric(usedata$feature)
  usedata <- usedata[which(usedata[,"treatment"] %in% c("E2","ED","EDTAM","EDFULV")),]
  
  E2dataWT <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  EDdataWT <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  TAMdataWT <- usedata[which(usedata[,"treatment"] == "EDTAM" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  FULVdataWT <- usedata[which(usedata[,"treatment"] == "EDFULV" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  
  E2dataYS <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
  EDdataYS <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
  TAMdataYS <- usedata[which(usedata[,"treatment"] == "EDTAM" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
  FULVdataYS <- usedata[which(usedata[,"treatment"] == "EDFULV" & usedata[,"celltype"] == "MCF7YS") ,"feature"]

  WT_E2ED  <- twosidePVttestDIR(E2dataWT,EDdataWT)
  WT_E2TAM <- twosidePVttestDIR(E2dataWT,TAMdataWT)
  WT_E2FULV <- twosidePVttestDIR(E2dataWT,FULVdataWT)
  WT_EDTAM  <- twosidePVttestDIR(EDdataWT,TAMdataWT)
  WT_EDFULV <- twosidePVttestDIR(EDdataWT,FULVdataWT)
  YS_E2ED  <- twosidePVttestDIR(E2dataYS,EDdataYS)
  YS_E2TAM <- twosidePVttestDIR(E2dataYS,TAMdataYS)
  YS_E2FULV <- twosidePVttestDIR(E2dataYS,FULVdataYS)
  YS_EDTAM <- twosidePVttestDIR(EDdataYS,TAMdataYS)
  YS_EDFULV <- twosidePVttestDIR(EDdataYS,FULVdataYS)
  
  allPV <- rbind(allPV, c(thisFeature, WT_E2ED  ,WT_E2TAM ,WT_E2FULV,WT_EDTAM ,WT_EDFULV,YS_E2ED  ,YS_E2TAM ,YS_E2FULV,YS_EDTAM ,YS_EDFULV))
  df <- data.frame(Method = c(rep("ED", length(EDdataWT)), rep("E2", length(E2dataWT)),rep("TAM", length(TAMdataWT)),rep("FULV", length(FULVdataWT))),
                   Accuracy = c(EDdataWT, E2dataWT,TAMdataWT,FULVdataWT))
  desired_order <- c("ED", "E2", "TAM","FULV")
  df$Method <- factor(df$Method, levels = desired_order)
  colors <- rep("#A11EFF",4)#c("E2" = "red", "ED" = "blue")
  p1 <- ggplot(df, aes(x = Method, y = Accuracy, fill = Method)) +
    geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
    stat_summary(fun = median, geom = "crossbar", width = 0.25, color = "black", size = 0.5) +  # Add horizontal line for median
    scale_fill_manual(values = colors) +  # Custom fill colors for violins
    theme_minimal() +  guides(fill = "none") + # Remove color legend
    labs(title = "MCF7 WT", y = thisFeature, x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1))

  df <- data.frame(Method = c(rep("ED", length(EDdataYS)), rep("E2", length(E2dataYS)),rep("TAM", length(TAMdataYS)),rep("FULV", length(FULVdataYS))),
                   Accuracy = c(EDdataYS, E2dataYS,TAMdataYS,FULVdataYS))
  desired_order <- c("ED", "E2", "TAM","FULV")
  df$Method <- factor(df$Method, levels = desired_order)
  colors <- rep("#1026FF",4)#c("E2" = "red", "ED" = "blue")
  p2 <- ggplot(df, aes(x = Method, y = Accuracy, fill = Method)) +
    geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
    stat_summary(fun = median, geom = "crossbar", width = 0.25, color = "black", size = 0.5) +  # Add horizontal line for median
    scale_fill_manual(values = colors) +  # Custom fill colors for violins
    theme_minimal() +  guides(fill = "none") + # Remove color legend
    labs(title = "MCF7 YS", y = thisFeature, x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1))
  pdf(file=paste0("cbBatch_featureBox/violin_",thisFeature,".pdf"),width=8,height=6)
  grid.arrange(p1,p2, ncol=2)
  dev.off()
  
  pdf(file=paste0("cbBatch_featureBox/boxplot_",thisFeature,".pdf"),width=8,height=6)
  par(mfrow=c(1,2),mar=c(4,4,2,2))
  boxplot(EDdataWT,E2dataWT,TAMdataWT,FULVdataWT,
          names = c("ED","E2","TAM","FULV"), las=2, col="#A11EFF",
          outline = F,ylab=thisFeature, main="MCF7 WT")
  boxplot(EDdataYS,E2dataYS,TAMdataYS,FULVdataYS,
          names = c("ED","E2","TAM","FULV"), las=2, col="#1026FF",
          outline = F,ylab=thisFeature, main="MCF7 YS")
  dev.off()
}


# ---------------------------------
# Feature Analysis (batch-separated)
# ---------------------------------
#### batch-separated analysis, boxplot comparing features between E2 and ED
GGvioplot <- function(E2data, EDdata,outname,YL){
  data1 <- E2data
  data2 <- EDdata
  df <- data.frame(Method = c(rep("E2", length(data1)), rep("ED", length(data2))),
                   Accuracy = c(data1, data2))
  desired_order <- c("E2", "ED")
  df$Method <- factor(df$Method, levels = desired_order)
  colors <- c("E2" = "red", "ED" = "blue")
  
  p <- ggplot(df, aes(x = Method, y = Accuracy, fill = Method)) +
    geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
    stat_summary(fun = median, geom = "crossbar", width = 0.25, color = "black", size = 0.5) +  # Add horizontal line for median
    scale_fill_manual(values = colors) +  # Custom fill colors for violins
    theme_minimal() +
    labs(title = outname, y = YL, x = "Methods") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1))
  ggsave(paste0("sepBatchE2ED_featureBox/violin_MCF7WT_sepBatchE2ED_",YL,"_",outname,".pdf"))
  pdf(file=paste0("sepBatchE2ED_featureBox/boxplot_MCF7WT_sepBatchE2ED_",YL,"_",outname,".pdf"))
  boxplot(E2data, EDdata,outline=F,col=c("red","blue"),names=c("E2","ED"),ylab=YL,main=outname,xlab=twosidePVttestEX(E2data,EDdata))
  dev.off()
}
colnames(allPV) <- c("cmp","WT_E2ED","WT_E2TAM","WT_E2FULV","WT_EDTAM","WT_EDFULV","YS_E2ED","YS_E2TAM","YS_E2FULV","YS_EDTAM","YS_EDFULV")
write.table(allPV,file="allPV.txt",row.names=F,col.names=T,sep="\t",quote=F)

for(thisFeature in useFeatures){
  for(databatch in sort(unique(meta_mat_use[,1]))[1:6]){
    useidx <- which(meta_mat_use[,1] == databatch)
    usedata<- as.data.frame( cbind(meta_mat_use[useidx,], domain_use[useidx,thisFeature]) )
    colnames(usedata) <- c("name","celltype","chan1","chan2","batch","treatment","rep","cell","feature")
    usedata$feature <- as.numeric(usedata$feature)
    E2data <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    EDdata <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    GGvioplot(E2data,EDdata,databatch,thisFeature)
  }
}

#### batch-separated analysis, barplot for T-statistics score
#### material for Figure S1C
outT <- c()
outPV <- c()
for(thisFeature in useFeatures){
  for(databatch in sort(unique(meta_mat_use[,1]))){
    useidx <- which(meta_mat_use[,1] == databatch)
    usedata<- as.data.frame( cbind(meta_mat_use[useidx,], domain_use[useidx,thisFeature]) )
    colnames(usedata) <- c("name","celltype","chan1","chan2","batch","treatment","rep","cell","feature")
    usedata$feature <- as.numeric(usedata$feature)
    E2data <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    EDdata <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    EDTAMdata<- usedata[which(usedata[,"treatment"] == "EDTAM" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    EDFULVdata <- usedata[which(usedata[,"treatment"] == "EDFULV" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
    E2dataYS <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
    EDdataYS <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
    EDTAMdataYS <- usedata[which(usedata[,"treatment"] == "EDTAM" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
    EDFULVdataYS <- usedata[which(usedata[,"treatment"] == "EDFULV" & usedata[,"celltype"] == "MCF7YS") ,"feature"]
    
    WT_E2ED <- twosidePVttestEX(E2data,EDdata)
    WT_EDTAMED <- twosidePVttestEX(EDTAMdata,EDdata)
    WT_EDFULVED <- twosidePVttestEX(EDFULVdata,EDdata)
    WT_E2EDTAM <- twosidePVttestEX(E2data,EDTAMdata)
    YS_E2ED <- twosidePVttestEX(E2dataYS,EDdataYS)
    YS_EDTAMED <- twosidePVttestEX(EDTAMdataYS,EDdataYS)
    YS_EDFULVED <- twosidePVttestEX(EDFULVdataYS,EDdataYS)
    YS_E2EDTAM <- twosidePVttestEX(E2dataYS,EDTAMdataYS)
    outPV <- rbind(outPV, c(thisFeature, databatch, WT_E2ED[1], WT_EDTAMED[1], WT_EDFULVED[1], WT_E2EDTAM[1],YS_E2ED[1], YS_EDTAMED[1], YS_EDFULVED[1], YS_E2EDTAM[1]))
    outT <- rbind(outT, c(thisFeature, databatch, WT_E2ED[2], WT_EDTAMED[2], WT_EDFULVED[2], WT_E2EDTAM[2],YS_E2ED[2], YS_EDTAMED[2], YS_EDFULVED[2], YS_E2EDTAM[2]))
  }
}
colnames(outPV) <- c("feature","batch","WT_E2ED","WT_EDTAMED","WT_EDFULVED","WT_E2EDTAM","YS_E2ED","YS_EDTAMED","YS_EDFULVED","YS_E2EDTAM")
write.table(outPV, file=paste0("domainLevel_sepBatch_PV.txt"),row.names=F,col.names=T,sep="\t",quote=F)
colnames(outT) <- c("feature","batch","WT_E2ED","WT_EDTAMED","WT_EDFULVED","WT_E2EDTAM","YS_E2ED","YS_EDTAMED","YS_EDFULVED","YS_E2EDTAM")
write.table(outT, file=paste0("domainLevel_sepBatch_Tscore.txt"),row.names=F,col.names=T,sep="\t",quote=F)


# ploting
Tdata <- read.table("domainLevel_sepBatch_Tscore.txt",header=T)
PVdata <-read.table("domainLevel_sepBatch_PV.txt",header=T)

barbox <- function(values, p_values, XL, M){
  #print(t.test(values,mu=0)$p.val)
  #print(-log10(t.test(values,mu=0)$p.val))
  Tvs0 <- round(-log10(t.test(values,mu=0)$p.val),2)
  Tvs0_label <- as.character(Tvs0)
  df <- data.frame(Index = factor(1:length(values)), Values = values, p_values = p_values)
  df$significance <- ifelse(df$p_values < 0.05, "Significant", "Not Significant")
  plotobj <- ggplot(df, aes(x = Index, y = Values)) +
    geom_bar(stat = "identity", aes(fill = significance), width = 0.5) +
    geom_boxplot(aes(x = as.factor(length(values) + 1), y = Values), width = 0.5, fill = "grey") +
    geom_jitter(aes(x = as.factor(length(values) + 1), y = Values, color = significance), width = 0.2, size = 3) +
    scale_x_discrete(limits = c(as.character(1:length(values)), as.character(length(values) + 1)),labels = c(XL, Tvs0_label)) +
    scale_fill_manual(values = c("Significant" = "darkblue", "Not Significant" = "lightblue")) +
    scale_color_manual(values = c("Significant" = "darkblue", "Not Significant" = "lightblue")) +
    labs(x = "", y = "t-statistics", title = M) +
    theme_minimal() +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "none")
  
  return(plotobj) 
}

# Sphericity
plot1 <- barbox(Tdata[1:6,"WT_E2ED"],     PVdata[1:6,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "Sphericity WT E2vsED")
plot2 <- barbox(Tdata[1:6,"WT_EDTAMED"],  PVdata[1:6,"WT_EDTAMED"],  c("B2","B3","B4","B5","B6","B7"), "Sphericity WT EDTAMvsED")
plot3 <- barbox(Tdata[1:6,"WT_EDFULVED"], PVdata[1:6,"WT_EDFULVED"], c("B2","B3","B4","B5","B6","B7"), "Sphericity WT EDFULVvsED")
plot4 <- barbox(Tdata[1:6,"WT_E2EDTAM"],  PVdata[1:6,"WT_E2EDTAM"],  c("B2","B3","B4","B5","B6","B7"), "Sphericity WT E2vsEDTAM")
plot5 <- barbox(Tdata[7:13,"YS_E2ED"],    PVdata[7:13,"YS_E2ED"],    c("B1","B2","B3","B4","B5","B6","B7"), "Sphericity YS E2vsED")
plot6 <- barbox(Tdata[7:13,"YS_EDTAMED"], PVdata[7:13,"YS_EDTAMED"], c("B1","B2","B3","B4","B5","B6","B7"), "Sphericity YS EDTAMvsED")
plot7 <- barbox(Tdata[7:13,"YS_EDFULVED"],PVdata[7:13,"YS_EDFULVED"],c("B1","B2","B3","B4","B5","B6","B7"), "Sphericity YS EDFULVvsED")
plot8 <- barbox(Tdata[7:13,"YS_E2EDTAM"], PVdata[7:13,"YS_E2EDTAM"], c("B1","B2","B3","B4","B5","B6","B7"), "Sphericity YS E2vsEDTAM")
plots <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf("sepBatch_bar_Sphericity.pdf", width = 16, height = 8)  # Adjust the width and height as necessary
grid.arrange(grobs = plots, nrow = 2, ncol = 4)
dev.off()

# corH3K27acER
plot1 <- barbox(Tdata[14:19,"WT_E2ED"],     PVdata[14:19,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "corH3K27ac&ER WT E2vsED")
plot2 <- barbox(Tdata[14:19,"WT_EDTAMED"],  PVdata[14:19,"WT_EDTAMED"],  c("B2","B3","B4","B5","B6","B7"), "corH3K27ac&ER WT EDTAMvsED")
plot3 <- barbox(Tdata[14:19,"WT_EDFULVED"], PVdata[14:19,"WT_EDFULVED"], c("B2","B3","B4","B5","B6","B7"), "corH3K27ac&ER WT EDFULVvsED")
plot4 <- barbox(Tdata[14:19,"WT_E2EDTAM"],  PVdata[14:19,"WT_E2EDTAM"],  c("B2","B3","B4","B5","B6","B7"), "corH3K27ac&ER WT E2vsEDTAM")
plot5 <- barbox(Tdata[20:26,"YS_E2ED"],     PVdata[20:26,"YS_E2ED"],    c("B1","B2","B3","B4","B5","B6","B7"),  "corH3K27ac&ER YS E2vsED")
plot6 <- barbox(Tdata[20:26,"YS_EDTAMED"],  PVdata[20:26,"YS_EDTAMED"], c("B1","B2","B3","B4","B5","B6","B7"),  "corH3K27ac&ER YS EDTAMvsED")
plot7 <- barbox(Tdata[20:26,"YS_EDFULVED"], PVdata[20:26,"YS_EDFULVED"],c("B1","B2","B3","B4","B5","B6","B7"),  "corH3K27ac&ER YS EDFULVvsED")
plot8 <- barbox(Tdata[20:26,"YS_E2EDTAM"],  PVdata[20:26,"YS_E2EDTAM"], c("B1","B2","B3","B4","B5","B6","B7"),  "corH3K27ac&ER YS E2vsEDTAM")
plots <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf("sepBatch_bar_corH3K27acER.pdf", width = 16, height = 8)  # Adjust the width and height as necessary
grid.arrange(grobs = plots, nrow = 2, ncol = 4)
dev.off()

# volumn
plot1 <- barbox(Tdata[27:32,"WT_E2ED"],     PVdata[27:32,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "volumn WT E2vsED")
plot2 <- barbox(Tdata[27:32,"WT_EDTAMED"],  PVdata[27:32,"WT_EDTAMED"],  c("B2","B3","B4","B5","B6","B7"), "volumn WT EDTAMvsED")
plot3 <- barbox(Tdata[27:32,"WT_EDFULVED"], PVdata[27:32,"WT_EDFULVED"], c("B2","B3","B4","B5","B6","B7"), "volumn WT EDFULVvsED")
plot4 <- barbox(Tdata[27:32,"WT_E2EDTAM"],  PVdata[27:32,"WT_E2EDTAM"],  c("B2","B3","B4","B5","B6","B7"), "volumn WT E2vsEDTAM")
plot5 <- barbox(Tdata[33:39,"YS_E2ED"],     PVdata[33:39,"YS_E2ED"],    c("B1","B2","B3","B4","B5","B6","B7"),  "volumn YS E2vsED")
plot6 <- barbox(Tdata[33:39,"YS_EDTAMED"],  PVdata[33:39,"YS_EDTAMED"], c("B1","B2","B3","B4","B5","B6","B7"),  "volumn YS EDTAMvsED")
plot7 <- barbox(Tdata[33:39,"YS_EDFULVED"], PVdata[33:39,"YS_EDFULVED"],c("B1","B2","B3","B4","B5","B6","B7"),  "volumn YS EDFULVvsED")
plot8 <- barbox(Tdata[33:39,"YS_E2EDTAM"],  PVdata[33:39,"YS_E2EDTAM"], c("B1","B2","B3","B4","B5","B6","B7"),  "volumn YS E2vsEDTAM")
plots <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf("sepBatch_bar_volumn.pdf", width = 16, height = 8)  # Adjust the width and height as necessary
grid.arrange(grobs = plots, nrow = 2, ncol = 4)
dev.off()

# EXvolumn
plot1 <- barbox(Tdata[40:45,"WT_E2ED"],     PVdata[40:45,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "EXvolumn WT E2vsED")
plot2 <- barbox(Tdata[40:45,"WT_EDTAMED"],  PVdata[40:45,"WT_EDTAMED"],  c("B2","B3","B4","B5","B6","B7"), "EXvolumn WT EDTAMvsED")
plot3 <- barbox(Tdata[40:45,"WT_EDFULVED"], PVdata[40:45,"WT_EDFULVED"], c("B2","B3","B4","B5","B6","B7"), "EXvolumn WT EDFULVvsED")
plot4 <- barbox(Tdata[40:45,"WT_E2EDTAM"],  PVdata[40:45,"WT_E2EDTAM"],  c("B2","B3","B4","B5","B6","B7"), "EXvolumn WT E2vsEDTAM")
plot5 <- barbox(Tdata[46:52,"YS_E2ED"],     PVdata[46:52,"YS_E2ED"],    c("B1","B2","B3","B4","B5","B6","B7"),  "EXvolumn YS E2vsED")
plot6 <- barbox(Tdata[46:52,"YS_EDTAMED"],  PVdata[46:52,"YS_EDTAMED"], c("B1","B2","B3","B4","B5","B6","B7"),  "EXvolumn YS EDTAMvsED")
plot7 <- barbox(Tdata[46:52,"YS_EDFULVED"], PVdata[46:52,"YS_EDFULVED"],c("B1","B2","B3","B4","B5","B6","B7"),  "EXvolumn YS EDFULVvsED")
plot8 <- barbox(Tdata[46:52,"YS_E2EDTAM"],  PVdata[46:52,"YS_E2EDTAM"], c("B1","B2","B3","B4","B5","B6","B7"),  "EXvolumn YS E2vsEDTAM")
plots <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf("sepBatch_bar_EXvolumn.pdf", width = 16, height = 8)  # Adjust the width and height as necessary
grid.arrange(grobs = plots, nrow = 2, ncol = 4)
dev.off()

# dist boundary
plot1 <- barbox(Tdata[53:58,"WT_E2ED"],     PVdata[53:58,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "distBoundary WT E2vsED")
plot2 <- barbox(Tdata[53:58,"WT_EDTAMED"],  PVdata[53:58,"WT_EDTAMED"],  c("B2","B3","B4","B5","B6","B7"), "distBoundary WT EDTAMvsED")
plot3 <- barbox(Tdata[53:58,"WT_EDFULVED"], PVdata[53:58,"WT_EDFULVED"], c("B2","B3","B4","B5","B6","B7"), "distBoundary WT EDFULVvsED")
plot4 <- barbox(Tdata[53:58,"WT_E2EDTAM"],  PVdata[53:58,"WT_E2EDTAM"],  c("B2","B3","B4","B5","B6","B7"), "distBoundary WT E2vsEDTAM")
plot5 <- barbox(Tdata[59:65,"YS_E2ED"],     PVdata[59:65,"YS_E2ED"],    c("B1","B2","B3","B4","B5","B6","B7"),  "distBoundary YS E2vsED")
plot6 <- barbox(Tdata[59:65,"YS_EDTAMED"],  PVdata[59:65,"YS_EDTAMED"], c("B1","B2","B3","B4","B5","B6","B7"),  "distBoundary YS EDTAMvsED")
plot7 <- barbox(Tdata[59:65,"YS_EDFULVED"], PVdata[59:65,"YS_EDFULVED"],c("B1","B2","B3","B4","B5","B6","B7"),  "distBoundary YS EDFULVvsED")
plot8 <- barbox(Tdata[59:65,"YS_E2EDTAM"],  PVdata[59:65,"YS_E2EDTAM"], c("B1","B2","B3","B4","B5","B6","B7"),  "distBoundary YS E2vsEDTAM")
plots <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf("sepBatch_bar_distBoundary.pdf", width = 16, height = 8)  # Adjust the width and height as necessary
grid.arrange(grobs = plots, nrow = 2, ncol = 4)
dev.off()

#### summary of T-statistics scores
#### material for Figure 2B,D,F
boxonly <- function(values, p_values, XL, M, y1,y2){
  Tvs0 <- round(-log10(t.test(values,mu=0)$p.val),2)
  Tvs0_label <- as.character(Tvs0)
  df <- data.frame(Index = factor(1:length(values)), Values = values, p_values = p_values)
  df$significance <- ifelse(df$p_values < 0.05, "Significant", "Not Significant")
  plotobj <- ggplot(df, aes(x = Index, y = Values)) +
    geom_boxplot(aes(x = 1, y = Values), width = 0.2, fill = "grey",outlier.shape = NA) +
    geom_jitter(aes(x = 1, y = Values, color = significance), width = 0.2, size = 3) +
    scale_color_manual(values = c("Significant" = "lightblue", "Not Significant" = "lightblue")) +
    labs(x = "", y = "t-statistics", title = M) +
    ylim(y1, y2) +   xlim(0, 4)+
    theme_minimal() + geom_hline(yintercept = 0, color = "blue") +  # Horizontal line at y=0
    theme(legend.position = "none")
  return(plotobj) 
}

plot1 <- boxonly(Tdata[1:6,"WT_E2ED"],     PVdata[1:6,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "Sphericity WT E2vsED",-15,15)
plot2 <- boxonly(Tdata[14:19,"WT_E2ED"],     PVdata[14:19,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "corH3K27ac&ER WT E2vsED",-20,20)
plot3 <- boxonly(Tdata[27:32,"WT_E2ED"],     PVdata[27:32,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "volumn WT E2vsED",-45,45)
plot4 <- boxonly(Tdata[53:58,"WT_E2ED"],     PVdata[53:58,"WT_E2ED"],     c("B2","B3","B4","B5","B6","B7"), "distBoundary WT E2vsED",-21,21)
pdf("sepBatch_summaryBox.pdf", width = 8, height = 8)  # Adjust the width and height as necessary
plots <- list(plot1, plot2, plot3, plot4)
grid.arrange(grobs = plots, nrow = 2, ncol = 2)
dev.off()


# ---------------------------------
# Other backup analysis
# ---------------------------------
for(thisFeature in useFeatures ){
  usedata<- as.data.frame( cbind(meta_mat_use, domain_use[,thisFeature]) )
  colnames(usedata) <- c("name","celltype","chan1","chan2","batch","treatment","rep","cell","feature")
  usedata$feature <- as.numeric(usedata$feature)
  usedata <- usedata[which(usedata[,"treatment"] %in% c("E2","ED","EDTAM","EDFULV")),]
  
  E2data <- usedata[which(usedata[,"treatment"] == "E2" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  EDdata <- usedata[which(usedata[,"treatment"] == "ED" & usedata[,"celltype"] == "MCF7WT") ,"feature"]
  WT_E2ED <- twosidePVttestEX(E2data,EDdata)
  
  data1 <- E2data
  data2 <- EDdata
  df <- data.frame(Method = c(rep("E2", length(data1)), rep("ED", length(data2))),
                   Accuracy = c(data1, data2))
  desired_order <- c("E2", "ED")
  df$Method <- factor(df$Method, levels = desired_order)
  colors <- c("E2" = "red", "ED" = "blue")
  
  p <- ggplot(df, aes(x = Method, y = Accuracy, fill = Method)) +
    geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
    stat_summary(fun = median, geom = "crossbar", width = 0.25, color = "black", size = 0.5) +  # Add horizontal line for median
    scale_fill_manual(values = colors) +  # Custom fill colors for violins
    theme_minimal() +
    labs(title = "MCF7WT all batches", y = thisFeature, x = "Methods") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1))
  ggsave(paste0("violin_MCF7WT_cbBatchE2ED_",thisFeature,".pdf"))
  pdf(file=paste0("boxplot_MCF7WT_cbBatchE2ED_",thisFeature,".pdf"))
  #boxplot(E2data, EDdata,outline=F,col=c("red","blue"),names=c("E2","ED"),ylab=thisFeature,main="MCF7WT all batches")
  #  scale_y_continuous(limits = c(-30, 30)) +  # Setting the y limits here
  boxplot(E2data, EDdata,outline=F,col=c("red","blue"),names=c("E2","ED"),
          ylab=thisFeature,xlab=c(paste0("PV= ",WT_E2ED[1]),paste0("t= ",WT_E2ED[2])),main="MCF7WT all batches")
  dev.off()
}
