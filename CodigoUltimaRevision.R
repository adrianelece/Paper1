library(DSS)
library(tidyverse)
library(data.table)
library(ChIPseeker)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(clusterProfiler)
library (ggupset)
library (ggimage)
library(ggpubr)
library(org.Bt.eg.db)
library(GenomicRanges)
library(ReactomePA)
library(ChIPpeakAnno)
library(ggVennDiagram)
library(plotly)
library(ggridges)
library(peRReo)
library(wesanderson)
library(treemap)
library(ggvenn)

txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene
paleta = latin_palette("badbunny2",11,type = "continuous")
paleta1 = paleta[11:1]
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May")

chip4x = fread("4xRegsToAnno.txt")
chip10x = fread("10xRegsToAnno.txt")


##### 1x
peak = readPeakFile("1xRegsToAnno.txt")

peakAnno1x <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno1x) + 
  scale_fill_manual(values = latin_palette("badbunny2",11,type = "continuous"))
peakAnnodf1x=data.frame(peakAnno1x)

peakAnnodf_tss1x = peakAnnodf1x[which(abs(peakAnnodf1x$distanceToTSS)<=50000),]
Regsless50kb = paste0(gsub("chr","",peakAnnodf_tss1x$seqnames),":",peakAnnodf_tss1x$start)
write.table(Regsless50kb,"RegsLess50kb1x.txt",quote = F,col.names = F,row.names = F)


######## 4x
peak = readPeakFile("4xRegsToAnno.txt")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno) + 
  scale_fill_manual(values = latin_palette("badbunny2",11,type = "continuous"))
peakAnnodf=data.frame(peakAnno)

peakAnnodf_tss = peakAnnodf[which(abs(peakAnnodf$distanceToTSS)<=50000),]
Regsless50kb = paste0(gsub("chr","",peakAnnodf_tss$seqnames),":",peakAnnodf_tss$start)
write.table(Regsless50kb,"RegsLess50kb4x.txt",quote = F,col.names = F,row.names = F)

######## 7x
peak7 = readPeakFile("7xRegsToAnno.txt")

peakAnno7 <- annotatePeak(peak7, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno7) + 
  scale_fill_manual(values = latin_palette("badbunny2",11,type = "continuous"))
peakAnnodf7=data.frame(peakAnno7)

peakAnnodf_tss7 = peakAnnodf7[which(abs(peakAnnodf7$distanceToTSS)<=50000),]
Regsless50kb7 = paste0(gsub("chr","",peakAnnodf_tss7$seqnames),":",peakAnnodf_tss7$start)
write.table(Regsless50kb7,"RegsLess50kb7x.txt",quote = F,col.names = F,row.names = F)

##### 10x
peak10 = readPeakFile("10xRegsToAnno.txt")

peakAnno10 <- annotatePeak(peak10, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno10)
peakAnnodf10=data.frame(peakAnno10)

peakAnnodf_tss10 = peakAnnodf10[which(abs(peakAnnodf10$distanceToTSS)<=50000),]
Regsless50kb10 = paste0(gsub("chr","",peakAnnodf_tss10$seqnames),":",peakAnnodf_tss10$start)
write.table(Regsless50kb10,"RegsLess50kb10x.txt",quote = F,col.names = F,row.names = F)

##### Plots

peakAnno1x@annoStat %>% group_by(Feature) %>% arrange(Frequency) %>% mutate(key = as.factor(paste0(Feature," (",round(Frequency,2),")"))) %>% 
  ggplot(.,aes(x="", y=Frequency, fill=key)) +
  scale_fill_manual(values=paleta1)+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  ggtitle("CpGs called by RRBS at a minimum coverage of 1x") +
  theme(text = element_text(size = 25),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Genomic feature"))


peakAnno@annoStat %>% group_by(Feature) %>% arrange(Frequency) %>% mutate(key = as.factor(paste0(Feature," (",round(Frequency,2),")"))) %>% 
  ggplot(.,aes(x="", y=Frequency, fill=key)) +
  scale_fill_manual(values=paleta1)+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  ggtitle("CpGs called by RRBS at a minimum coverage of 4x") +
  theme(text = element_text(size = 25),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Genomic feature"))


peakAnno7@annoStat %>% group_by(Feature) %>% arrange(Frequency) %>% mutate(key = as.factor(paste0(Feature," (",round(Frequency,2),")"))) %>% 
  ggplot(.,aes(x="", y=Frequency, fill=key)) +
  scale_fill_manual(values=paleta1)+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  ggtitle("CpGs called by RRBS at a minimum coverage of 7x") +
  theme(text = element_text(size = 25),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Genomic feature"))

peakAnno10@annoStat %>% group_by(Feature) %>% arrange(Frequency) %>% mutate(key = as.factor(paste0(Feature," (",round(Frequency,2),")"))) %>% 
ggplot(.,aes(x="", y=Frequency, fill=key)) +
  scale_fill_manual(values=paleta1)+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  ggtitle("CpGs called by RRBS at a minimum coverage of 10x") +
  theme(text = element_text(size = 25),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Genomic feature"))



##### DM analyses


##################
#####  RRBS ######
##################


setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/10x/")
filestodo1x = list.files("")
files1x = filestodo1x[grep("DSS_clean*",filestodo1x)]

for(i in 1:length(files1x)){
  assign(paste0("meth_",gsub("_synthese_CpG_filtered.txt","",gsub("DSS_clean_SQM2_","",files1x[i]))),fread(files1x[i]))
}

colnames(meth_10_15) = c("chr", "pos", "N", "X")
colnames(meth_4_14  ) = c("chr", "pos", "N", "X")
colnames(meth_4_23  ) = c("chr", "pos", "N", "X")
colnames(meth_8_2   ) = c("chr", "pos", "N", "X")
colnames(meth_8_21) = c("chr", "pos", "N", "X")
colnames(meth_8_19  ) = c("chr", "pos", "N", "X")

BSobj1x = makeBSseqData(list(meth_10_15,
                             meth_4_14,
                             meth_4_23,
                             meth_8_2,
                             meth_8_21,
                             meth_8_19  
                             ),
                        c("10_15",
                          "4_14",
                          "4_23",
                          "8_2",
                          "8_21",
                          "8_19"
                        ))

BSobj1x

dmlTest.sm1x = DMLtest(BSobj1x, group1=c("4_14","4_23","8_21"),
                       group2=c("8_2","8_19","10_15"), 
                       smoothing=TRUE)
head(dmlTest.sm1x)

#write.table(dmlTest.sm1x,"DMLtest_1x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)
dmls1x = callDML(dmlTest.sm1x, p.threshold=0.05)
head(dmls1x)
write.table(dmls1x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls21x = callDML(dmlTest.sm1x, delta=0.2, p.threshold=0.05)
head(dmls21x)
#write.table(dmls21x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs1x = callDMR(dmlTest.sm1x, p.threshold=0.05)
head(dmrs1x)
write.table(dmrs1x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr1x = dmls1x[which(dmls1x$fdr<=0.05),]
#write.table(dmlsfdr1x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta1x = dmls21x[which(dmls21x$fdr<=0.05),]
#write.table(dmlsfdrdelta1x,"dml005delta2fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)



##### 4x

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x")
filestodo4x = list.files("")
files4x = filestodo4x[grep("DSS_clean*",filestodo4x)]

for(i in 1:length(files4x)){
  assign(paste0("meth_",gsub("_synthese_CpG_filtered.txt","",gsub("DSS_clean_SQM2_","",files4x[i]))),fread(files4x[i]))
}

colnames(meth_10_15) = c("chr", "pos", "N", "X")
colnames(meth_4_14  ) = c("chr", "pos", "N", "X")
colnames(meth_4_23  ) = c("chr", "pos", "N", "X")
colnames(meth_8_2   ) = c("chr", "pos", "N", "X")
colnames(meth_8_21) = c("chr", "pos", "N", "X")
colnames(meth_8_19  ) = c("chr", "pos", "N", "X")

BSobj4x = makeBSseqData(list(meth_10_15,
                             meth_4_14,
                             meth_4_23,
                             meth_8_2,
                             meth_8_21,
                             meth_8_19  
),
c("10_15",
  "4_14",
  "4_23",
  "8_2",
  "8_21",
  "8_19"
))

BSobj4x

dmlTest.sm4x = DMLtest(BSobj4x, group1=c("4_14","4_23","8_21"),
                       group2=c("8_2","8_19","10_15"), 
                       smoothing=TRUE)
head(dmlTest.sm4x)

write.table(dmlTest.sm4x,"DMLtest_4x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)
dmls4x = callDML(dmlTest.sm4x, p.threshold=0.05)
head(dmls4x)
write.table(dmls4x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls24x = callDML(dmlTest.sm4x, delta=0.2, p.threshold=0.05)
head(dmls24x)
write.table(dmls24x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs4x = callDMR(dmlTest.sm4x, p.threshold=0.05)
head(dmrs4x)
write.table(dmrs4x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr4x = dmls4x[which(dmls4x$fdr<=0.05),]
write.table(dmlsfdr4x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta4x = dmls24x[which(dmls24x$fdr<=0.05),]
write.table(dmlsfdrdelta4x,"dml005delta2fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)


##### 7x

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x")
filestodo7x = list.files("")
files7x = filestodo[grep("DSS_clean*",filestodo)]

for(i in 1:length(files7x)){
  assign(paste0("meth_",gsub("_synthese_CpG_filtered.txt","",gsub("DSS_clean_SQM2_","",files7x[i]))),fread(files7x[i]))
}

colnames(meth_10_15) = c("chr", "pos", "N", "X")
colnames(meth_4_14  ) = c("chr", "pos", "N", "X")
colnames(meth_4_23  ) = c("chr", "pos", "N", "X")
colnames(meth_8_2   ) = c("chr", "pos", "N", "X")
colnames(meth_8_21) = c("chr", "pos", "N", "X")
colnames(meth_8_19  ) = c("chr", "pos", "N", "X")

BSobj7x = makeBSseqData(list(meth_10_15,
                             meth_4_14,
                             meth_4_23,
                             meth_8_2,
                             meth_8_21,
                             meth_8_19  
),
c("10_15",
  "4_14",
  "4_23",
  "8_2",
  "8_21",
  "8_19"
))

BSobj7x

dmlTest.sm7x = DMLtest(BSobj7x, group1=c("4_14","4_23","8_21"),
                       group2=c("8_2","8_19","10_15"), 
                       smoothing=TRUE)
head(dmlTest.sm7x)

write.table(dmlTest.sm7x,"DMLtest_7x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)
dmls7x = callDML(dmlTest.sm7x, p.threshold=0.05)
head(dmls7x)
write.table(dmls7x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls27x = callDML(dmlTest.sm7x, delta=0.2, p.threshold=0.05)
head(dmls27x)
write.table(dmls27x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs7x = callDMR(dmlTest.sm7x, p.threshold=0.05)
head(dmrs7x)
write.table(dmrs7x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr7x = dmls7x[which(dmls7x$fdr<=0.05),]
write.table(dmlsfdr7x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta7x = dmls27x[which(dmls27x$fdr<=0.05),]
write.table(dmlsfdrdelta7x,"dml005delta2fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)



##### 10x

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/10x")
filestodo10x = list.files("")
files10x = filestodo[grep("DSS_clean*",filestodo)]

for(i in 1:length(files10x)){
  assign(paste0("meth_",gsub("_synthese_CpG_filtered.txt","",gsub("DSS_clean_SQM2_","",files10x[i]))),fread(files10x[i]))
}

colnames(meth_10_15) = c("chr", "pos", "N", "X")
colnames(meth_4_14  ) = c("chr", "pos", "N", "X")
colnames(meth_4_23  ) = c("chr", "pos", "N", "X")
colnames(meth_8_2   ) = c("chr", "pos", "N", "X")
colnames(meth_8_21) = c("chr", "pos", "N", "X")
colnames(meth_8_19  ) = c("chr", "pos", "N", "X")

BSobj10x = makeBSseqData(list(meth_10_15,
                              meth_4_14,
                              meth_4_23,
                              meth_8_2,
                              meth_8_21,
                              meth_8_19  
),
c("10_15",
  "4_14",
  "4_23",
  "8_2",
  "8_21",
  "8_19"
))

BSobj10x

dmlTest.sm10x = DMLtest(BSobj10x, group1=c("4_14","4_23","8_21"),
                        group2=c("8_2","8_19","10_15"), 
                        smoothing=TRUE)
head(dmlTest.sm10x)

write.table(dmlTest.sm10x,"DMLtest_10x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)
dmls10x = callDML(dmlTest.sm10x, p.threshold=0.05)
head(dmls10x)
write.table(dmls10x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls210x = callDML(dmlTest.sm10x, delta=0.2, p.threshold=0.05)
head(dmls210x)
write.table(dmls210x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs10x = callDMR(dmlTest.sm10x, p.threshold=0.05)
head(dmrs10x)
write.table(dmrs10x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr10x = dmls10x[which(dmls10x$fdr<=0.05),]
write.table(dmlsfdr10x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta10x = dmls210x[which(dmls210x$fdr<=0.05),]
write.table(dmlsfdrdelta10x,"dml005delta2fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)



######################
######################
###    Graficos    ###
######################
######################

###RRBS DMC todos

dmls4x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x/dml005.tsv")
dmls7x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x/dml005.tsv")


### RRBS1x

dmlsfdr1x = fread("1x/dml005fdr005.tsv")
dmlsfdrdelta1x = fread("1x/dml005delta2fdr005.tsv")

ggplot() +
  geom_density(aes(diff, fill = "DMC1x"), alpha = .2, data = dmls1x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC1xfdr"), alpha = .2, data = dmlsfdr1x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC1xdelta"), alpha = .2, data = dmlsfdrdelta1x,linetype="dotdash",color="blue")+
  scale_fill_manual(name = "DMC", values = c(DMC1x = "red", DMC1xfdr = "green", DMC1xdelta = "blue")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (1x)") + 
  ylab("Density")

### RRBS4x

dmlsfdr4x = fread("4x/dml005fdr005.tsv")
dmlsfdrdelta4x = fread("4x/dml005delta2fdr005.tsv")

ggplot() +
  geom_density(aes(diff, fill = "DMC4x"), alpha = .2, data = dmls4x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC4xfdr"), alpha = .2, data = dmlsfdr4x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC4xdelta"), alpha = .2, data = dmlsfdrdelta4x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMC", values = c(DMC4x = "red", DMC4xfdr = "green", DMC4xdelta = "blue")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (4x)") + 
  ylab("Density")

### RRBS7x
dmlsfdr7x = fread("7x/dml005fdr005.tsv")
dmlsfdrdelta7x = fread("7x/dml005delta2fdr005.tsv")


ggplot() +
  geom_density(aes(diff, fill = "DMC7x", color = "DMC7x"), alpha = .2, data = dmls7x, linetype = "dashed") +
  geom_density(aes(diff, fill = "DMC7xfdr", color = "DMC7xfdr"), alpha = .2, data = dmlsfdr7x, linetype = "dotted") +
  geom_density(aes(diff, fill = "DMC7xdelta", color = "DMC7xdelta"), alpha = .2, data = dmlsfdrdelta7x, linetype = "dotdash") +
  scale_fill_manual(name = "DMC", values = c(DMC7x = "red", DMC7xfdr = "green", DMC7xdelta = "blue")) +
  scale_color_manual(name = "DMC", values = c(DMC7x = "red", DMC7xfdr = "green", DMC7xdelta = "blue")) +
  guides(fill = guide_legend(), color = guide_legend()) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (7x)") +
  ylab("Density")


### RRBS10x
ggplot() +
  geom_density(aes(diff, fill = "DMC10x"), alpha = .2, data = dmls10x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC10xfdr"), alpha = .2, data = dmlsfdr10x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC10xdelta"), alpha = .2, data = dmlsfdrdelta10x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMC", values = c(DMC10x = "red", DMC10xfdr = "green", DMC10xdelta = "blue")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (10x)") + 
  ylab("Density")

### RRBS todos

ggplot() +
  geom_density(aes(diff, fill = "DMC1x", color = "DMC1x"), alpha = .2, data = dmlsfdrdelta1x, linetype = "dashed") +
  geom_density(aes(diff, fill = "DMC4x", color = "DMC4x"), alpha = .2, data = dmlsfdrdelta4x, linetype = "dotted") +
  geom_density(aes(diff, fill = "DMC7x", color = "DMC7x"), alpha = .2, data = dmlsfdrdelta7x, linetype = "dotdash") +
  scale_fill_manual(name = "DMC", values = c(DMC1x = "red", DMC4x = "green", DMC7x = "blue")) +
  scale_color_manual(name = "DMC", values = c(DMC1x = "red", DMC4x = "green", DMC7x = "blue")) +
  guides(fill = guide_legend(), color = guide_legend()) +
  theme_minimal() +
  xlab("Difference between the estimated means. Filtered by FDR and δ=2") +
  ylab("Density")


ggplot() +
  geom_density(aes(diff.Methy, fill = "DMR1x"), alpha = .2, data = dmrs1x,linetype="dashed",color="red") +
  geom_density(aes(diff.Methy, fill= "DMR4x"), alpha = .2, data = dmrs4x,linetype="dotted",color="green") +
  geom_density(aes(diff.Methy, fill= "DMR7x"), alpha = .2, data = dmrs7x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMR", values = c(DMR1x = "red", DMR4x = "green", DMR7x = "blue")) +
  theme_minimal() +
  xlab("Difference between the estimated means of DMRs") + 
  ylab("Density")

### ONT COMPARACION

ontdmls4x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/4x/dml005.tsv")
ontdmls7x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/7x/dml005.tsv")
cols = wes_palette("Chevalier1",4)
breaks =seq(-1, 1, by=0.25)
labels = as.character(breaks)

comp1 = ggplot() +
  geom_density(aes(diff, fill = "RRBS DMC4x", linetype = "RRBS DMC4x"), alpha = 0.5, data = dmls4x, color = cols[1]) +
  geom_density(aes(diff, fill = "RRBS DMC7x", linetype = "RRBS DMC7x"), alpha = 0.5, data = dmls7x, color = cols[2]) +
  geom_density(aes(diff, fill = "ONT DMC4x", linetype = "ONT DMC4x"), alpha = 1, data = ontdmls4x, color = cols[3]) +
  geom_density(aes(diff, fill = "ONT DMC7x", linetype = "ONT DMC7x"), alpha = 0.5, data = ontdmls7x, color = cols[4]) +
  scale_fill_manual(name = "DMC", values = c("RRBS DMC4x" = cols[1], "RRBS DMC7x" = cols[2], "ONT DMC4x" = cols[3], "ONT DMC7x" = cols[4])) +
  scale_linetype_manual(name = "DMC", values = c("RRBS DMC4x" = "solid", "RRBS DMC7x" = "solid", "ONT DMC4x" = "solid", "ONT DMC7x" = "solid")) +
  geom_vline(xintercept=-0.2, linetype="dashed",color="red",size=0.75)+
  geom_vline(xintercept=0.2, linetype="dashed",color="red",size=0.75)+
  theme_classic2() +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key = element_rect(fill = "black"),
    plot.title = element_text(hjust = 0.5)) +
  xlab("Difference between the estimated means") +
  ylab("Density") +
  ggtitle("DMC filtered by p-value")+
  scale_x_continuous(limits = c(-1, 1), breaks = breaks, labels = labels)


# RRBS dmlsfdrdelta4x
# RRBS dmlsfdrdelta7x

ontdmlsfdrdelta4x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/4x/dml005delta2fdr005.tsv")
ontdmlsfdrdelta7x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/7x/dml005delta2fdr005.tsv")
cols = wes_palette("Chevalier1",4)


comp2 = ggplot() +
  geom_density(aes(diff, fill = "RRBS DMC FDR and Delta4x", linetype = "RRBS DMC FDR and Delta4x"), alpha = 0.5, data = dmlsfdrdelta4x, color = cols[1]) +
  geom_density(aes(diff, fill = "RRBS DMC FDR and Delta7x", linetype = "RRBS DMC FDR and Delta7x"), alpha = 0.5, data = dmlsfdrdelta7x, color = cols[2]) +
  geom_density(aes(diff, fill = "ONT DMC FDR and Delta4x", linetype = "ONT DMC FDR and Delta4x"), alpha = 1, data = ontdmlsfdrdelta4x, color = cols[3]) +
  geom_density(aes(diff, fill = "ONT DMC FDR and Delta7x", linetype = "ONT DMC FDR and Delta7x"), alpha = 0.5, data = ontdmlsfdrdelta7x, color = cols[4]) +
  scale_fill_manual(name = "DMC", values = c("RRBS DMC FDR and Delta4x" = cols[1], "RRBS DMC FDR and Delta7x" = cols[2], "ONT DMC FDR and Delta4x" = cols[3], "ONT DMC FDR and Delta7x" = cols[4])) +
  scale_linetype_manual(name = "DMC", values = c("RRBS DMC FDR and Delta4x" = "solid", "RRBS DMC FDR and Delta7x" = "solid", "ONT DMC FDR and Delta4x" = "solid", "ONT DMC FDR and Delta7x" = "solid")) +
  geom_vline(xintercept=-0.2, linetype="dashed",color="red",size=0.75)+
  geom_vline(xintercept=0.2, linetype="dashed",color="red",size=0.75)+
  theme_classic2() +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key = element_rect(fill = "black"),
    plot.title = element_text(hjust = 0.5)) +
    xlab("Difference between the estimated means") +
    ylab("Density") +
  ggtitle("\n DMC filtered by p-value and FDR of 0.05 and delta of 0.2")+
scale_x_continuous(limits = c(-1, 1), breaks = breaks, labels = labels)



allcomps = ggarrange(comp1 + rremove("xlab") + rremove("ylab"),
                     comp2 + rremove("xlab") + rremove("ylab"),ncol=1,
                     common.legend = TRUE, legend="bottom"
                   )

allcompsanno = annotate_figure(allcomps, 
                              left = textGrob("Density\n", rot = 90, vjust = 1, gp = gpar(cex = 2))
                              )



#########################
######  Annotation  #####
#########################

###########4x###########

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x")
chip4x = data.frame(dmlsfdrdelta4x$chr,dmlsfdrdelta4x$pos,dmlsfdrdelta4x$pos+1,paste0("r",1:nrow(dmlsfdrdelta4x)))
write.table(chip4x,"Chip4x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip4x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)
peakAnnodftss = peakAnnodf[which(abs(peakAnnodf$distanceToTSS)<=50000),]
genes4xrrbs = unique(peakAnnodftss$SYMBOL)
write.table(genes4x,"GenesRRBS4x.tsv",sep="\t",quote=F,row.names = F,col.names = F)

############
##1x####
#########
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/1x")

chip1x = data.frame(dmlsfdrdelta1x$chr,dmlsfdrdelta1x$pos,dmlsfdrdelta1x$pos+1,paste0("r",1:nrow(dmlsfdrdelta1x)))
write.table(chip1x,"Chip1x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip1x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

genes1x = unique(peakAnnodf_tss$SYMBOL)

############
##7x####
#########
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x")

chip7x = data.frame(dmlsfdrdelta7x$chr,dmlsfdrdelta7x$pos,dmlsfdrdelta7x$pos+1,paste0("r",1:nrow(dmlsfdrdelta7x)))
write.table(chip7x,"Chip7x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip7x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)
peakAnnodftss = peakAnnodf[which(abs(peakAnnodf$distanceToTSS)<=50000),]
genes7xrrbs = unique(peakAnnodftss$SYMBOL)
write.table(genes7x,"GenesRRBS7x.tsv",sep="\t",quote=F,row.names = F,col.names = F)

#library("gprofiler2")
#gostres <- gost(query = genes7x, 
#                organism = "btaurus", ordered_query = FALSE, 
#                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
#                measure_underrepresentation = FALSE, evcodes = FALSE, 
#                user_threshold = 0.05, correction_method = "g_SCS", 
#                domain_scope = "annotated", custom_bg = NULL, 
#                numeric_ns = "", sources = NULL, as_short_link = FALSE)


############
##10x####
#########
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/10x")

chip10x = data.frame(dmlsfdrdelta10x$chr,dmlsfdrdelta10x$pos,dmlsfdrdelta10x$pos+1,paste0("r",1:nrow(dmlsfdrdelta10x)))
write.table(chip10x,"Chip10x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip10x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

genes10x = unique(peakAnnodf_tss$SYMBOL)


genes = list(RRBS4x=genes4x,RRBS7x=genes7x)

ggvenn(
  genes, 
  fill_color = wes_palette("Chevalier1")[1:4],
  stroke_size = 0.5, set_name_size = 5)





######### ONT ############

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo")

ontdmlsfdr4x = fread("4x/dml005delta2fdr005.tsv")
ontdmlsfdr7x = fread("7x/dml005delta2fdr005.tsv")

cols=wes_palette("Chevalier1",4)


ggplot() +
  geom_density(aes(diff, fill = "RRBS4x"), alpha = .2, data = dmlsfdrdelta1x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "RRBS7x"), alpha = .2, data = dmlsfdrdelta4x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "ONT4x"), alpha = .2, data = ontdmlsfdr4x,linetype="dotdash",color="blue") +
  geom_density(aes(diff, fill= "ONT7x"), alpha = .2, data = ontdmlsfdr7x,linetype="dotdash",color="pink") +
  scale_fill_manual(name = "DMC", values = c(RRBS4x = cols[1], RRBS7x = cols[2], ONT4x = cols[3],ONT7x = cols[4])) +
  theme_minimal() +
  xlab("Difference between the estimated means. Filtered by FDR and δ=2") + 
  ylab("Density")



ggplot() +
  geom_density(aes(diff, fill = "RRBS4x", color = "RRBS4x"), alpha = .5, data = dmlsfdrdelta1x, linetype = "dashed") +
  geom_density(aes(diff, fill = "RRBS7x", color = "RRBS7x"), alpha = .5, data = dmlsfdrdelta4x, linetype = "dotted") +
  geom_density(aes(diff, fill = "ONT4x", color = "ONT4x"), alpha = .8, data = ontdmlsfdr4x, linetype = "dotdash") +
  geom_density(aes(diff, fill = "ONT7x", color = "ONT7x"), alpha = .8, data = ontdmlsfdr7x, linetype = "dotdash") +
  scale_fill_manual(name = "DMC", values = c(RRBS4x = cols[1], RRBS7x = cols[2], ONT4x = cols[3],ONT7x = cols[4])) +
  scale_color_manual(name = "DMC", values = c(RRBS4x = cols[1], RRBS7x = cols[2], ONT4x = cols[3],ONT7x = cols[4])) +
  guides(fill = guide_legend(), color = guide_legend()) +
  theme_minimal() +
  xlab("Difference between the estimated means. Filtered by FDR and δ=2") +
  ylab("Density")

chip4xont = data.frame(paste0("chr",ontdmlsfdr4x$chr),ontdmlsfdr4x$pos,ontdmlsfdr4x$pos+1,paste0("r",1:nrow(ontdmlsfdr4x)))
write.table(chip4x,"Chip4xONT.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip4xONT.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

peakAnnodftss = peakAnnodf[which(abs(peakAnnodf$distanceToTSS)<=50000),]
genes4x = unique(peakAnnodftss$SYMBOL)
write.table(genes4x,"GenesONT4x.tsv",sep="\t",quote=F,row.names = F,col.names = F)

chip7xont = data.frame(paste0("chr",ontdmlsfdr7x$chr),ontdmlsfdr7x$pos,ontdmlsfdr7x$pos+1,paste0("r",1:nrow(ontdmlsfdr7x)))
write.table(chip7xont,"Chip7xONT.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip7xONT.tsv")

peakAnno7x <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno7x)
peakAnnodf7x=data.frame(peakAnno7x)
peakAnnodftss = peakAnnodf[which(abs(peakAnnodf7x$distanceToTSS)<=50000),]

genes7x = unique(peakAnnodftss$SYMBOL)
write.table(genes7x,"GenesONT7x.tsv",sep="\t",quote=F,row.names = F,col.names = F)



### Anno comp

files = list(RRBS4x = "C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x/Chip4x.tsv",
             RRBS7x = "C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x/Chip7x.tsv",
             ONT4x  = "C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/4x/Chip4x.tsv",
             ONT7x  = "C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/7x/Chip7x.tsv"
             )

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
tagHeatmap(tagMatrixList,xlim=c(-3000, 3000))

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)


names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
vennplot(genes)

plotAnnoBar(peakAnnoList)+ggtitle("Genetic features associated to the DMCs") + 
  scale_fill_manual(values=latin_palette("badbunny2", n = 10, type = "continuous")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position = "bottom")
plotDistToTSS(peakAnnoList)+ggtitle("Distance from the DMCs to the TSS")+ 
  scale_fill_manual(values=latin_palette("badbunny2", n = 6, type = "continuous")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position = "bottom")


#genes = list(ONT4x=genes4x,ONT7x=genes7x,RRBS4x=genes4xrrbs,RRBS7x=genes7xrrbs)
ggVennDiagram(genes)+
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")

ggvenn(
  genes, 
  fill_color = wes_palette("Chevalier1")[1:4],
  stroke_size = 0.5, set_name_size = 5
)

genes$RRBS4x = genes$RRBS4x[!duplicated(genes$RRBS4x)]
genes$RRBS7x = genes$RRBS7x[!duplicated(genes$RRBS7x)]
genes$ONT4x = genes$ONT4x[!duplicated(genes$ONT4x)]
genes$ONT7x = genes$ONT7x[!duplicated(genes$ONT7x)]

library(eulerr)
plot(euler(genes, shape = "ellipse"),fills=wes_palette("Chevalier1")[1:4], labels = c("RRBS 4x","RRBS 7x","ONT 4x","ONT 7x"), quantities = TRUE)




########### QTLS

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/qtls")

rrbsqtl4x = fread("QTLsIn_trim_4xRRBS.tsv")
rrbsqtl4x[,1:4]$V4 = gsub(".[(].*[)]", "", rrbsqtl4x[,1:4]$V4)
rrbsqtl4x[,1:4]$V4 = gsub("QTL","",rrbsqtl4x[,1:4]$V4)
rrbsqtl4x[,1:4]$V4 = trimws(rrbsqtl4x[,1:4]$V4)

rrbsqtl7x = fread("QTLsIn_trim_7xRRBS.tsv")
rrbsqtl7x[,1:4]$V4 = gsub(".[(].*[)]", "", rrbsqtl7x[,1:4]$V4)
rrbsqtl7x[,1:4]$V4 = gsub("QTL","",rrbsqtl7x[,1:4]$V4)
rrbsqtl7x[,1:4]$V4 = trimws(rrbsqtl7x[,1:4]$V4)

ontqtls4x = fread("QTLsIn_trim_4xONT.tsv")
ontqtls4x = ontqtls4x[,1:4]
ontqtls4x$V4 = gsub(".[(].*[)]", "", ontqtls4x$V4)
ontqtls4x$V4 = gsub("QTL","",ontqtls4x$V4)
ontqtls4x$V4 = trimws(ontqtls4x$V4)

ontqtls7x = fread("QTLsIn_trim_7xONT.tsv")
ontqtls7x = ontqtls7x[,1:4]
ontqtls7x$V4 = gsub(".[(].*[)]", "", ontqtls7x$V4)
ontqtls7x$V4 = gsub("QTL","",ontqtls7x$V4)
ontqtls7x$V4 = trimws(ontqtls7x$V4)

qtlsbed = fread("QTL_trimmed.bed")
qtlsbed = qtlsbed[,1:4]
qtlsbed$V4 = gsub(".[(].*[)]", "", qtlsbed$V4)
qtlsbed$V4 = gsub("QTL","",qtlsbed$V4)
qtlsbed$V4 = trimws(qtlsbed$V4)


catrrbs4x = as.data.frame(table(rrbsqtl4x$V4))
catrrbs4x$tech = "RRBS 4x"

catrrbs7x = as.data.frame(table(rrbsqtl7x$V4))
catrrbs7x$tech = "RRBS 7x"

catontqtls4x = as.data.frame(table(ontqtls4x$V4))
catontqtls4x$tech = "ONT 4x"

catontqtls7x = as.data.frame(table(ontqtls7x$V4))
catontqtls7x$tech = "ONT 7x"

catqtlsbed = as.data.frame(table(qtlsbed$V4))
catqtlsbed$tech = "DataBase"

qtlscats = rbind(catrrbs4x,catrrbs7x,catontqtls4x,catontqtls7x,catqtlsbed)
qtlscats_trim = qtlscats[-which(qtlscats$Freq<7),]

data_new <- qtlscats_trim                              # Replicate data
data_new$tech <- factor(data_new$tech,      # Reordering group factor levels
                        levels = c("RRBS 4x", "RRBS 7x", "ONT 4x", "ONT 7x", "DataBase"))


cols = wes_palette("Chevalier1", 4)

data_new %>% filter(tech !="DataBase") %>%
ggplot(., aes(fill=tech, y=reorder(Var1,-Freq), x=Freq)) +
  geom_bar(position="dodge", stat="identity",col="black") +
  scale_fill_manual(values=cols) +
  ggtitle("") +
  facet_wrap(~tech,ncol=4) +
  theme_bw()+
  theme(legend.position="none",
        strip.background =element_rect(colour="black", fill="white"),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 45, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 25, margin = margin())
  ) +
  xlab("\nNumber of DMCs associated with the QTL")+
  ylab("QTL")



data_new %>% filter(tech == "DataBase") %>% filter(Freq >800) %>%
  ggplot(., aes(fill=tech, y=reorder(Var1,-Freq), x=Freq)) +
  geom_bar(position="dodge", stat="identity",col="white") +
  scale_fill_manual(values="black") +
  ggtitle("") +
  facet_wrap(~tech,ncol=4) +
  theme_bw()+
  theme(legend.position="none",
        strip.background =element_rect(colour="black", fill="white"),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 45, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 25, margin = margin())
  ) +
  xlab("\nNumber of DMCs associated with the QTL")+
  ylab("QTL")


####KDE

peakRRBS4x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x/Chip4x.tsv")
peakRRBS7x = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x/Chip7x.tsv")
peakONT4x  = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/4x/Chip4x.tsv")
peakONT7x  = fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/7x/Chip7x.tsv")


colnames(peakRRBS4x)=c("chr","start","end","rx")
colnames(peakRRBS7x)=c("chr","start","end","rx")
colnames(peakONT4x)=c("chr","start","end","rx")
colnames(peakONT7x)=c("chr","start","end","rx")


peakrrbs2 = peakRRBS4x %>% group_by(chr) %>% mutate(Lenght = max(end)-min(start), Breaks=round(Lenght/60))
peakrrbs2$window = NA

for(i in 1:nrow(peakrrbs2)){
  for(j in 1:60){
    if(peakrrbs2$start[i]<=peakrrbs2$Breaks[i]*j&&is.na(peakrrbs2$window[i])){
      peakrrbs2$window[i] = j
    }
  }
}
peakrrbs2$window[which(is.na(peakrrbs2$window))]=60


peakrrbs3 = peakRRBS7x %>% group_by(chr) %>% mutate(Lenght = max(end)-min(start), Breaks=round(Lenght/60))
peakrrbs3$window = NA

for(i in 1:nrow(peakrrbs3)){
  for(j in 1:60){
    if(peakrrbs3$start[i]<=peakrrbs3$Breaks[i]*j&&is.na(peakrrbs3$window[i])){
      peakrrbs3$window[i] = j
    }
  }
}
peakrrbs3$window[which(is.na(peakrrbs3$window))]=60



peakont2 = peakONT4x %>% group_by(chr) %>% mutate(Lenght = max(end)-min(start), Breaks=round(Lenght/60))
peakont2$window = NA

for(i in 1:nrow(peakont2)){
  for(j in 1:60){
    if(peakont2$start[i]<=peakont2$Breaks[i]*j&&is.na(peakont2$window[i])){
      peakont2$window[i] = j
    }
  }
}
peakont2$window[which(is.na(peakont2$window))]=60


peakont3 = peakONT7x %>% group_by(chr) %>% mutate(Lenght = max(end)-min(start), Breaks=round(Lenght/60))
peakont3$window = NA

for(i in 1:nrow(peakont3)){
  for(j in 1:60){
    if(peakont3$start[i]<=peakont3$Breaks[i]*j&&is.na(peakont3$window[i])){
      peakont3$window[i] = j
    }
  }
}
peakont3$window[which(is.na(peakont3$window))]=60



peakrrbs2$tech = "RRBS4x"
peakrrbs3$tech = "RRBS7x"
peakont2$tech = "ONT4x"
peakont3$tech = "ONT7x"

peakstot = rbind(peakrrbs2,peakrrbs3,peakont2,peakont3)

library(wesanderson)
peakstot %>%
  mutate(tech = factor(tech, levels = c("RRBS4x","RRBS7x","ONT4x","ONT7x"))) %>%
  ggplot(aes(x = window, fill = tech)) +
  geom_density(position = 'identity', alpha = 0.8, adjust = 1) +
  guides(fill=guide_legend(title="Technique")) +
  scale_fill_manual(values = wes_palette("Chevalier1")) +
  facet_wrap(~tech)+
  ylab("Density of the marks \n") + 
  xlab("\n Relative position of the marks along the chromosome") +
  scale_x_continuous(breaks=c(2, 31, 59), 
                     labels=c('Chr start', '50% chr length', 'Chr end')) +
  theme_bw() +
  theme(text = element_text(size=25),
        legend.position = "bottom",
        strip.background=element_rect(colour="black",
                                      fill="white"),
        axis.text.x = element_text(size=15))

##### Enrichment
compGO<- compareCluster(geneCluster   = genes,
                        OrgDb = org.Bt.eg.db,
                        fun           = "enrichGO",
                        ont = "ALL",
                        pool = T,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH")

dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")




######### HEATMAP
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/heatmap")
heatm = fread("Heatmap.tsv")

heatfil = heatm %>% filter(Coverage!="10x")


corrs1 = heatfil[,c(1,3,4,5)]
corrs2 = heatfil[,c(1,3,6,7)]
corrs1$status = "Filtered"
corrs2$status = "Unfiltered"
corrs3 = rbind(corrs1,corrs2,use.names=F)
corrs3$comstatus = paste0(corrs3$Coverage,"_",corrs3$status)
colnames(corrs3)=c("Sample","Coverage","Correlation","Sites","Filt","Status")



corplot = corrs3 %>%
  mutate(Status = factor(Status, levels = c("10x_Unfiltered","7x_Unfiltered","4x_Unfiltered",
                                            "10x_Filtered","7x_Filtered","4x_Filtered"
  ))) %>%
  mutate(Sample = fct_reorder(Sample,desc(Sample))) %>% 
  ggplot(aes(x=Sample, y = Correlation, fill=Status)) +
  geom_bar(position="dodge", stat="identity", alpha=1, width=.8)+
  coord_flip() +
  xlab("") +
  ylim(0,1)+
  scale_fill_manual("Status", 
                    values = c("4x_Filtered"    = latin_palette("badbunny1")[1],
                               "7x_Filtered"    = latin_palette("badbunny1")[2],
                               "10x_Filtered"   = latin_palette("badbunny1")[3],
                               "4x_Unfiltered"  = latin_palette("badbunny1")[4],
                               "7x_Unfiltered"  = latin_palette("badbunny1")[5],
                               "10x_Unfiltered" = latin_palette("badbunny1")[7])) + 
  theme_bw() + 
  theme(text = element_text(size=20),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  facet_wrap(~Filt, scales = "free")

sitesplot = corrs3 %>%
  mutate(Status = factor(Status, levels = c("7x_Unfiltered","4x_Unfiltered",
                                            "7x_Filtered","4x_Filtered"
  ))) %>%
  mutate(Sample = fct_reorder(Sample,desc(Sample))) %>% 
  mutate(Sites = as.numeric(Sites)) %>%
  ggplot(aes(x=Sample, y = Sites, fill=Status)) +
  geom_bar(position="dodge", stat="identity", alpha=1, width=.8)+
  coord_flip() +
  xlab("") +
  scale_fill_manual("Status", 
                    values = c("4x_Filtered"    = latin_palette("badbunny1")[1],
                               "7x_Filtered"    = latin_palette("badbunny1")[2],
                               "4x_Unfiltered"  = latin_palette("badbunny1")[4],
                               "7x_Unfiltered"  = latin_palette("badbunny1")[5]
                               )) +
  theme_bw() + 
  theme(text = element_text(size=20),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  facet_wrap(~Filt, scales = "free")


segment_size <- 0.75

corrs3  %>% 
  ggplot(aes(y = Sample, x = Correlation)) +
  geom_point(col = "black") +
  geom_segment(aes(xend = 0, yend = Sample), col = "black", size = segment_size) +
  geom_text(
    aes(x = 0, label = Sample),
    size = 10,
    hjust = if_else(corrs3$Correlation > 0, 1, 0)
  ) +
  scale_y_discrete(breaks = NULL) +
  labs(x = 'Correlation', y = element_blank()) +
  theme_minimal() +
  facet_wrap(~Status)
  theme(
    panel.grid = element_line( linetype = 2, color = "red")
  )

  
########## HEATMAP
palden <- wes_palette("Zissou1", 15, type = "continuous")
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/ONTcovs/ONTcovs")
  

#### 4_14
s414 = fread("renamed_4_14_freq.tsv_4x.tsv")
colnames(s414) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x414 <- s414 %>% filter(covONT>=4&
                                  covBS>=4&
                                  freqBS>=0.1 & freqBS < 0.9&
                                  freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x414 <- s414 %>% filter(covONT>=7&
                                covBS>=7&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth414total = rbind(Meths4x414,Meths7x414)

corsmeth414 <- Meth414total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="4_14",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)


p414=Meth414total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
  #scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 1") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")


#### 4_23
s423 = fread("renamed_4_23_freq.tsv_4x.tsv")
colnames(s423) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x423 <- s423 %>% filter(covONT>=4&
                                covBS>=4&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x423 <- s423 %>% filter(covONT>=7&
                                covBS>=7&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth423total = rbind(Meths4x423,Meths7x423)

corsmeth423 <- Meth423total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="4_23",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)



p423=Meth423total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
  #scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 2") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")

#### 8_2
s82 = fread("renamed_8_2_freq.tsv_4x.tsv")
colnames(s82) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x82 <- s82 %>% filter(covONT>=4&
                              covBS>=4&
                              freqBS>=0.1 & freqBS < 0.9&
                              freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x82 <- s82 %>% filter(covONT>=7&
                              covBS>=7&
                              freqBS>=0.1 & freqBS < 0.9&
                              freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth82total = rbind(Meths4x82,Meths7x82)

corsmeth82 <- Meth82total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="8_2",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)



p82=Meth82total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
  #scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 3") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")

#### 8_19
s819 = fread("renamed_8_19_freq.tsv_4x.tsv")
colnames(s819) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x819 <- s819 %>% filter(covONT>=4&
                                covBS>=4&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x819 <- s819 %>% filter(covONT>=7&
                                covBS>=7&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth819total = rbind(Meths4x819,Meths7x819)

corsmeth819 <- Meth819total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="8_19",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)



p819=Meth819total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
  #scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 4") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")

#### 8_21
s821 = fread("renamed_8_21_2_freq.tsv_4x.tsv")
colnames(s821) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x821 <- s821 %>% filter(covONT>=4&
                                covBS>=4&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x821 <- s821 %>% filter(covONT>=7&
                                covBS>=7&
                                freqBS>=0.1 & freqBS < 0.9&
                                freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth821total = rbind(Meths4x821,Meths7x821)

corsmeth821 <- Meth821total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="8_21",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)



p821=Meth821total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
 # scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 5") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")

#### 10_15
s1015 = fread("renamed_10_15_freq.tsv_4x.tsv")
colnames(s1015) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
Meths4x1015 <- s1015 %>% filter(covONT>=4&
                                  covBS>=4&
                                  freqBS>=0.1 & freqBS < 0.9&
                                  freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "4x")

Meths7x1015 <- s1015 %>% filter(covONT>=7&
                                  covBS>=7&
                                  freqBS>=0.1 & freqBS < 0.9&
                                  freqONT>=0.1 & freqONT < 0.9) %>%
  mutate(Depth = "7x")

Meth1015total = rbind(Meths4x1015,Meths7x1015)

corsmeth1015 <- Meth1015total %>%
  group_by(Depth) %>% 
  mutate(Corrl = cor(freqBS,freqONT),Sample="10_15",Sites = n()) %>% 
  distinct(Corrl,Sample,Sites)



p1015=Meth1015total %>%
  mutate(across(Depth, factor, levels=c("4x","7x"))) %>%
  ggplot(aes(x=freqBS, y=freqONT)) +
  geom_density_2d_filled(alpha = 0.75,show.legend = NA) +
 # scale_fill_manual(values=palden)+
  facet_wrap(~Depth)+
  theme_bw() +
  xlab("Bisulfite Methylation Frequency") +
  ylab("ONT Methylation Frequency") +
  ggtitle("Sample 6") +
  theme(plot.title = element_text(hjust = 0.5 ),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        text = element_text(size=20),
        legend.position = "none")

allheat = ggarrange(p414 + rremove("ylab") + rremove("xlab"),
                    p423 + rremove("ylab") + rremove("xlab"),
                    p82 + rremove("ylab") + rremove("xlab"),
                    p819 + rremove("ylab") + rremove("xlab"),
                    p821 + rremove("ylab") + rremove("xlab"),
                    p1015 + rremove("ylab") + rremove("xlab"),
                    ncol = 2,nrow = 3,common.legend = F)

allheatanno = annotate_figure(allheat, 
                              left = textGrob("ONT methylation frequency\n", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                              bottom = textGrob("Bisulfite methylation frequency\n", gp = gpar(cex = 2)))



#####Comparison of sites
setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/CalledSites/")
sits = fread("SitesAndFilter.tsv")

sitslong = sits %>% 
  pivot_longer(`ONT 1x`:`Ratio 10x`, names_to = "Tech", values_to = "count")

remv = grep("Ratio",sitslong$Tech)

sitslongc = sitslong[-remv,]

sitsclean = sitslongc %>% group_by(Tech) %>% summarise(Mean=mean(count),SD=sd(count))

remv2 = grep("10x",sitsclean$Tech)
sitscleanc = sitsclean[-remv2,]
sitscleanc2 = sitscleanc
sitscleanc2$Tech <- factor(sitscleanc2$Tech,                                    # Change ordering manually
                  levels = c("ONT 1x", "RRBS 1x", "ONT 4x", "RRBS 4x", "ONT 7x", "RRBS 7x"))

cols = c(wes_palette("Rushmore1")[5],wes_palette("Rushmore1")[4],wes_palette("Chevalier1",4)[3],wes_palette("Chevalier1",4)[1],wes_palette("Chevalier1",4)[4],wes_palette("Chevalier1",4)[2])
cols = c(wes_palette("Rushmore1")[4],wes_palette("Rushmore1")[4],wes_palette("Chevalier1",4)[1],wes_palette("Chevalier1",4)[1],wes_palette("Chevalier1",4)[2],wes_palette("Chevalier1",4)[2])

ggplot(sitscleanc2, aes(x=Tech, y=Mean, fill=Tech)) + 
  geom_bar(color="black",stat="identity", alpha=1) +
  scale_fill_manual(values=cols) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2, colour="orange", alpha=0.9, size=1.3) +
  theme_bw() +
  xlab("") +
  ylab("Average number of CpGs called\n") +
  theme(text = element_text(size = 20),
        legend.position = "none")





### Anno repetido

rrbs4x=fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/4x/dml005delta2fdr005.tsv")
rrbs7x=fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/AnalisisPapers/Paper 1/May/DSS/7x/dml005delta2fdr005.tsv")
ont4x=fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/4x/dml005delta2fdr005.tsv")
ont7x=fread("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Estancia/Documents/semenmayo/marzo/marzo/7x/dml005delta2fdr005.tsv")

chip4xrrbs = data.frame(rrbs4x$chr,rrbs4x$pos,rrbs4x$pos+1,paste0("r",1:nrow(rrbs4x)))
chip7xrrbs = data.frame(rrbs7x$chr,rrbs7x$pos,rrbs7x$pos+1,paste0("r",1:nrow(rrbs7x)))
chip4xont = data.frame(paste0("chr",ont4x$chr),ont4x$pos,ont4x$pos+1,paste0("r",1:nrow(ont4x)))
chip7xont = data.frame(paste0("chr",ont7x$chr),ont7x$pos,ont7x$pos+1,paste0("r",1:nrow(ont7x)))

setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Escritorio/Papers/Paper1")
write.table(chip4xrrbs,"chiprrbs4x.tsv",sep="\t",quote=F,row.names = F,col.names = T)
write.table(chip7xrrbs,"chiprrbs7x.tsv",sep="\t",quote=F,row.names = F,col.names = T)
write.table(chip4xont,"chipont4x.tsv",sep="\t",quote=F,row.names = F,col.names = T)
write.table(chip7xont,"chipont7x.tsv",sep="\t",quote=F,row.names = F,col.names = T)
#rrbs4
peakr4x = readPeakFile("chiprrbs4x.tsv")
peakAnnorrbs4x <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Bt.eg.db")

peakAnnodfr4x=data.frame(peakAnnorrbs4x)
peakAnnodftssr4x = peakAnnodfr4x[which(abs(peakAnnodfr4x$distanceToTSS)<=50000),]
genes4xrrbs = unique(peakAnnodftssr4x$SYMBOL)

#rrbs7
peak = readPeakFile("chiprrbs7x.tsv")
peakAnno7xr <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Bt.eg.db")

peakAnnodfr7x=data.frame(peakAnno7xr)
peakAnnodftssr7x = peakAnnodfr7x[which(abs(peakAnnodfr7x$distanceToTSS)<=50000),]
genes7xrrbs = unique(peakAnnodftssr7x$SYMBOL)

#ont4
peako4 = readPeakFile("chipont4x.tsv")
peakAnnoo4 <- annotatePeak(peako4, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Bt.eg.db")

peakAnnodfo4=data.frame(peakAnnoo4)
peakAnnodftsso4 = peakAnnodfo4[which(abs(peakAnnodfo4$distanceToTSS)<=50000),]
genes4xont = unique(peakAnnodftsso4$SYMBOL)

#ont7
peako7 = readPeakFile("chipont7x.tsv")
peakAnno7xo <- annotatePeak(peako7, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Bt.eg.db")

peakAnnodfo7=data.frame(peakAnno7xo)
peakAnnodftsso7 = peakAnnodfo7[which(abs(peakAnnodfo7$distanceToTSS)<=50000),]
genes7xont = unique(peakAnnodftsso7$SYMBOL)


genes4xont=fread("C:/Users/alopez.catalina/Downloads/OneDrive_1_6-6-2023/genesONT4x.txt")
genes7xont=fread("C:/Users/alopez.catalina/Downloads/OneDrive_1_6-6-2023/genesONT7x.txt")


genes = list(ONT4x=genes4xont$X,ONT7x=genes7xont$X,RRBS4x=genes4xrrbs,RRBS7x=genes7xrrbs)
plot(euler(genes, shape = "ellipse"),fills=wes_palette("Chevalier1")[c(3,4,1,2)], quantities = TRUE)

ggVennDiagram(genes)+
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")

ggvenn(
  genes, 
  fill_color = wes_palette("Chevalier1")[1:4],
  stroke_size = 0.5, set_name_size = 5
)

test = fread("SampleXCFSFCUFSUF.tsv")

test %>% group_by(X) %>% summarise(mF=mean(SF),sdF=sd(SF))


library("gprofiler2")
gostres <- gost(query = genes7xont$X, 
                organism = "btaurus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)





