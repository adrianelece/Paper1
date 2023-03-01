# Paper1
## Preparacion de los datos para el DSS
Para evitar sobrepenalizacion hay que 
* Filtrar las regiones/citosinas que esten en el 100% de las muestras
* Quitar las que esten a mas de 50kb del TSS
* Pasar el DSS

## Filtrado
Para filtrar por el coverage de interes:
```
for file in *tsv; do echo $file; awk '($5>=XXX){print $0}' $file > "XXX_"$file; done
```
Donde XXX es el minimo coverage permitido.

## Crear regiones
Para poder seleccionar las citosinas en todas las muestras hay que crear regiones
```
for file in XXX_*; do echo $file; awk '{print $1":"$2"-"$3"\t"$7}' $file > 'network_'$file; done
```

## Crear lista con las regiones unicas
```
for file in network_*; do awk '{print $1}' $file > 'regions_'$file; done
cat regions_* > allregions.txt
sort allregions.txt | uniq > allregions.txt2
mv allregions.txt2 allregions.txt
```

## Extraer las regiones de la lista
```
python DatabyDir.py -r allregions.txt -s /folder/ -p "network*.tsv" -o AllSamples.tsv
```

## Clean the matrix
El output del programa anterior debe tener este aspecto 

```
10:100000318-100000330  1.0     .       1.0     .       .       1.0
        10:10000038-10000038    .       0.9     .       .       .       .
        10:100000448-100000462  1.0     0.75    1.0     1.0     .       1.0
        10:100001141-100001146  0.75    1.0     1.0     .       .       .
        10:100001344-100001354  0.667   1.0     1.0     1.0     .       .
        10:100002015-100002025  .       1.0     .       .       .       1.0
        10:10000309-10000309    0.667   1.0     .       .       .       .
        10:10000320-10000320    .       1.0     .       .       .       .
        10:100003781-100003781  1.0     .       .       .       .       .

```

Hay que limpiar las regiones con que tengan un . en alguna muestra. Para ello:

```
awk '($2!="."&&$3!="."&&$4!="."&&$5!="."&&$6!="."&&$7!="."){print $0}' AllSamples.tsv > RegsInAllSamples.txt
```

Y tiene que verse de la siguiente manera

```
        10:100007957-100007964  1.0     0.75    1.0     1.0     1.0     1.0
        10:100015958-100015971  0.0     0.0     0.167   1.0     0.0     0.333
        10:100017129-100017152  0.0     0.0     0.0     0.0     0.0     0.0
        10:100017181-100017208  0.0     0.0     0.0     0.0     0.0     0.0
        10:100024589-100024604  1.0     1.0     1.0     1.0     1.0     1.0
        10:10002736-10002744    0.0     0.0     0.0     0.0     0.0     0.0
        10:10002960-10002977    0.0     0.0     0.0     0.0     0.0     0.0
        10:10003027-10003037    0.0     0.0     0.0     0.0     0.0     0.0
        10:10003048-10003061    0.0     0.0     0.0     0.0     0.0     0.0

```

## Extraer las regiones para anotar

Estas regiones las anotamos para eliminar las que esten a mas de 50kb de los TSS

```
tail -n +2 RegsInAllSamples.txt | awk '{print $1}' | sed 's/:/\t/g' | sed 's/-/\t/g' | awk '{print "chr"$1"\t"$2-1"\t"$3+1"\t""r"NR}' > RegsToAnno.txt
```

## Annotation with ChipSeeker
```
library(DSS)
library(tidyverse)
library(data.table)
library(missForest)
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

txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene

setwd("XXXXX/")
peak = readPeakFile("RegsToAnno7x.txt")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

peakAnnodf_tss = peakAnnodf[which(peakAnnodf$distanceToTSS<=50000),]
Regsless50kb = paste0(gsub("chr","",peakAnnodf_tss$seqnames),":",peakAnnodf_tss$start)
write.table(Regsless50kb,"RegsLess50kb7x.txt",quote = F,col.names = F,row.names = F)

```


## Crear regiones en el archivo de methylation

```
for file in XXX_*; do echo $file; awk '{print $0"\t"$1":"$2"-"$3}' $file > "clean_"$file; done
```

## Crear el archivo para DSS
```
dos2unix RegsLess50kb7x.txt
for file in clean_*; do echo $file; grep -Fwf RegsLess50kb7x.txt $file | awk '{print "chr"$1"\t"$2"\t"$5"\t"$6}' > DSS/"DSS_"$file ; done
```

## DSS
```
library(DSS)
library(tidyverse)
library(data.table)
library(missForest)
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

txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene


setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo")
filestodo = list.files()
files = filestodo[grep("DSS_chip*",filestodo)]

for(i in 1:length(files)){
  assign(paste0("meth_",gsub(".tsvrenamed.tsv","",gsub("DSS_chip_clean_Region_4x_","",files[i]))),fread(files[i]))
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

#write.table(dmlTest.sm4x,"DMLtest_4x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls4x = callDML(dmlTest.sm4x, p.threshold=0.05)
head(dmls4x)
#write.table(dmls4x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls24x = callDML(dmlTest.sm4x, delta=0.2, p.threshold=0.05)
head(dmls24x)
#write.table(dmls24x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs4x = callDMR(dmlTest.sm4x, p.threshold=0.05)
head(dmrs4x)
#write.table(dmrs4x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr4x = dmls4x[which(dmls4x$fdr<=0.05),]
#write.table(dmlsfdr4x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta4x = dmls24x[which(dmls24x$fdr<=0.05),]
#write.table(dmlsfdrdelta4x,"dml005delta3fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)


#################
########7x#######
#################

setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/7x/")
filestodo = list.files()
files = filestodo[grep("*7x*",filestodo)]

for(i in 1:length(files)){
  assign(paste0("meth_",gsub(".tsvrenamed.tsv","",gsub("DSS_clean_7x_","",files[i]))),fread(files[i]))
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
write.table(dmlsfdrdelta7x,"dml005delta3fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)



###########
###10x#####
###########
setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/10x/")
filestodo = list.files()
files = filestodo[grep("*10*.tsv",filestodo)]

for(i in 1:length(files)){
  assign(paste0("meth_",gsub(".tsvrenamed.tsv","",gsub("DSS_clean_10_","",files[i]))),fread(files[i]))
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

#write.table(dmlTest.sm10x,"DMLtest_10x_smoothing.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls10x = callDML(dmlTest.sm10x, p.threshold=0.05)
head(dmls10x)
#write.table(dmls10x,"dml005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmls210x = callDML(dmlTest.sm10x, delta=0.2, p.threshold=0.05)
head(dmls210x)
#write.table(dmls210x,"dml005delta2.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmrs10x = callDMR(dmlTest.sm10x, p.threshold=0.05)
head(dmrs10x)
#write.table(dmrs10x,"dmrs005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdr10x = dmls10x[which(dmls10x$fdr<=0.05),]
#write.table(dmlsfdr10x,"dml005fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)

dmlsfdrdelta10x = dmls210x[which(dmls210x$fdr<=0.05),]
#write.table(dmls210x,"dml005delta3fdr005.tsv", sep="\t",quote = F, row.names = F,col.names = T)



######################
######################
###Probando cositas###
######################
######################
setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/imagenes")

tiff("DMC.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC4x"), alpha = .2, data = dmls4x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC7x"), alpha = .2, data = dmls7x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC10x"), alpha = .2, data = dmls10x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMC", values = c(DMC4x = "red", DMC7x = "green", DMC10x = "blue")) +
  theme_minimal() +
  theme(axis.text=element_text(size=12))+
  xlab("Difference between the estimated means") + 
  ylab("Density")
dev.off()

#ggplot(dmls4x, aes(x=diff, y=chr)) +
#  geom_density_ridges() +
#  theme_ridges() + 
#  theme(legend.position = "none")
tiff("Comp4x.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC4x"), alpha = .2, data = dmls4x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC4xfdr"), alpha = .2, data = dmlsfdr4x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC4xdelta"), alpha = .2, data = dmlsfdrdelta4x,linetype="dotdash",color="blue") +
  geom_density(aes(diff.Methy, fill= "DMR4x"), alpha = .2, data = dmrs4x,linetype="twodash",color="yellow") +
  scale_fill_manual(name = "DMC", values = c(DMC4x = "red", DMC4xfdr = "green", DMC4xdelta = "blue",DMR4x="yellow")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (4x)") + 
  ylab("Density")
dev.off()

tiff("Comp7x.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC7x"), alpha = .2, data = dmls7x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC7xfdr"), alpha = .2, data = dmlsfdr7x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC7xdelta"), alpha = .2, data = dmlsfdrdelta7x,linetype="dotdash",color="blue") +
  geom_density(aes(diff.Methy, fill= "DMR7x"), alpha = .2, data = dmrs7x,linetype="twodash",color="yellow") +
  scale_fill_manual(name = "DMC", values = c(DMC7x = "red", DMC7xfdr = "green", DMC7xdelta = "blue",DMR7x="yellow")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (7x)") + 
  ylab("Density")
dev.off()

tiff("Comp10x.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC10x"), alpha = .2, data = dmls10x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC10xfdr"), alpha = .2, data = dmlsfdr10x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC10xdelta"), alpha = .2, data = dmlsfdrdelta10x,linetype="dotdash",color="blue") +
  geom_density(aes(diff.Methy, fill= "DMR10x"), alpha = .2, data = dmrs10x,linetype="twodash",color="yellow") +
  scale_fill_manual(name = "DMC", values = c(DMC10x = "red", DMC10xfdr = "green", DMC10xdelta = "blue",DMR10x="yellow")) +
  theme_minimal() +
  xlab("Difference between the difference in the estimated means depending on the filter applied (10x)") + 
  ylab("Density")
dev.off()

tiff("Compfdr.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC4x"), alpha = .2, data = dmlsfdr4x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC7x"), alpha = .2, data = dmlsfdr7x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC10x"), alpha = .2, data = dmlsfdr10x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMC", values = c(DMC4x = "red", DMC7x = "green", DMC10x = "blue")) +
  theme_minimal() +
  xlab("Difference between the estimated. Filtered by FDR") + 
  ylab("Density")
dev.off()

tiff("Compfdrdelta02.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff, fill = "DMC4x"), alpha = .2, data = dmlsfdrdelta4x,linetype="dashed",color="red") +
  geom_density(aes(diff, fill= "DMC7x"), alpha = .2, data = dmlsfdrdelta7x,linetype="dotted",color="green") +
  geom_density(aes(diff, fill= "DMC10x"), alpha = .2, data = dmlsfdrdelta10x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMC", values = c(DMC4x = "red", DMC7x = "green", DMC10x = "blue")) +
  theme_minimal() +
  xlab("Difference between the estimated means. Filtered by FDR and δ=2") + 
  ylab("Density")
dev.off()

tiff("CompDMR.tiff", compression = "lzw",units = "px", width = 800, height = 589)
ggplot() +
  geom_density(aes(diff.Methy, fill = "DMR4x"), alpha = .2, data = dmrs4x,linetype="dashed",color="red") +
  geom_density(aes(diff.Methy, fill= "DMR7x"), alpha = .2, data = dmrs7x,linetype="dotted",color="green") +
  geom_density(aes(diff.Methy, fill= "DMR10x"), alpha = .2, data = dmrs10x,linetype="dotdash",color="blue") +
  scale_fill_manual(name = "DMR", values = c(DMR4x = "red", DMR7x = "green", DMR10x = "blue")) +
  theme_minimal() +
  xlab("Difference between the estimated means of DMRs") + 
  ylab("Density")
dev.off()



##################
####Annotation####
##################

###########4x###########
setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/")

chip4x = data.frame(dmlsfdrdelta4x$chr,dmlsfdrdelta4x$pos,dmlsfdrdelta4x$pos+1,paste0("r",1:nrow(dmlsfdrdelta4x)))
write.table(chip4x,"Chip4x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip4x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

genes4x = unique(peakAnnodf_tss$SYMBOL)


############
##7x####
#########
setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/7x/")

chip7x = data.frame(dmlsfdrdelta7x$chr,dmlsfdrdelta7x$pos,dmlsfdrdelta7x$pos+1,paste0("r",1:nrow(dmlsfdrdelta7x)))
write.table(chip7x,"Chip7x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip7x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

genes7x = unique(peakAnnodf_tss$SYMBOL)



############
##10x####
#########
setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/10x/")

chip10x = data.frame(dmlsfdrdelta10x$chr,dmlsfdrdelta10x$pos,dmlsfdrdelta10x$pos+1,paste0("r",1:nrow(dmlsfdrdelta10x)))
write.table(chip10x,"Chip10x.tsv",quote = F,col.names = F,row.names = F, sep = "\t")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peak = readPeakFile("Chip10x.tsv")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnno)
peakAnnodf=data.frame(peakAnno)

genes10x = unique(peakAnnodf_tss$SYMBOL)

##################
##Comp############
##################

setwd("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/comparacion")
rrbs = fread("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/enero/semensamples/Valentin/DSS-Diff_fertile_Subfertile-qvalue0.05.txt")
dssrrbb = data.frame(paste0("chr",rrbs$Chromosome),rrbs$Start-1,rrbs$End,paste0("r",1:nrow(rrbs)))
write.table(dssrrbb,"ChipSeekerRRBS.tsv",sep = "\t",quote = F, col.names = T,row.names = F)


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peakrrbs = readPeakFile("ChipSeekerRRBS.tsv")
peakont = readPeakFile("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/Chip4x.tsv")
peakont7x = readPeakFile("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/7x/Chip7x.tsv")
peakont10x = readPeakFile("C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/10x/Chip10x.tsv")

tiff("covRRBS.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
covplot(peakrrbs)+ggtitle("Distribution of the DMC detected by RRBS")
dev.off()
tiff("covONT4x.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
covplot(peakont)+ggtitle("Distribution of the DMC detected by ONT (4x)")
dev.off()
tiff("covONT7x.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
covplot(peakont7x)+ggtitle("Distribution of the DMC detected by ONT (7x)")
dev.off()
tiff("covONT10x.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
covplot(peakont10x)+ggtitle("Distribution of the DMC detected by ONT (7x)")
dev.off()

files=list(RRBS="ChipSeekerRRBS.tsv",ONT4x="C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/Chip4x.tsv",ONT7X="C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/7x/Chip7x.tsv",ONT10x="C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/semensamples/marzo/10x/Chip10x.tsv")
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
#plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
#tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
tiff("genomicregs.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
plotPeakProf2(files, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)
dev.off()
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

tiff("AnnoBar.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
plotAnnoBar(peakAnnoList)+ggtitle("Genetic features associated to the DMCs")
dev.off()

tiff("DistTTS.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
plotDistToTSS(peakAnnoList)+ggtitle("Distance from the DMCs to the TSS")
dev.off()


peakAnnorrbbs <- annotatePeak(peakrrbs, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnnorrbbs)
peakAnnodfrrbs=data.frame(peakAnnorrbbs)

peakAnnodf_tss_rrbs = peakAnnodfrrbs[which(peakAnnodfrrbs$distanceToTSS<=50000),]

genesrrbs = unique(peakAnnodf_tss_rrbs$SYMBOL)

peakAnnoont <- annotatePeak(peakont, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnnoont)
peakAnnodfont=data.frame(peakAnnoont)

peakAnnodf_tss_ont = peakAnnodfont[which(peakAnnodfont$distanceToTSS<=50000),]
genesont = unique(peakAnnodf_tss_ont$SYMBOL)

peakAnnoont7x <- annotatePeak(peakont7x, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnnoont7x)
peakAnnodfont7x=data.frame(peakAnnoont7x)

peakAnnodf_tss_ont7x = peakAnnodfont7x[which(peakAnnodfont7x$distanceToTSS<=50000),]
genesont7x = unique(peakAnnodf_tss_ont7x$SYMBOL)


peakAnnoont10x <- annotatePeak(peakont10x, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Bt.eg.db")
plotAnnoPie(peakAnnoont10x)
peakAnnodfont10x=data.frame(peakAnnoont10x)

peakAnnodf_tss_ont10x = peakAnnodfont10x[which(peakAnnodfont10x$distanceToTSS<=50000),]
genesont10x = unique(peakAnnodf_tss_ont10x$SYMBOL)

upsetplot(peakAnnoont, vennpie=T)
################Overlap de los genes unique ont 4x vs 7x vs rrbs
genes = list(ont=genesont,rrbs=genesrrbs,ont7x=genesont7x,ont10x=genesont10x)

tiff("venn2.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
ggVennDiagram(genes)+
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")
dev.off()

library(ggvenn)
library(wesanderson)
tiff("venn.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
ggvenn(
  genes, 
  fill_color = wes_palette("Chevalier1")[1:4],
  stroke_size = 0.5, set_name_size = 5
)
dev.off()

library(eulerr)
tiff("euler.tiff", width = 40, height = 20, pointsize = 1/300, units = 'cm', res = 300)
plot(euler(genes, shape = "ellipse"),fills=wes_palette("Chevalier1")[1:4], labels = c("ONT 4x","RRBS","ONT 7x","ONT 10x"), quantities = TRUE)
dev.off()
```


## QTL
QTL file for ARS-1.2 cattle genome can be downloaded from 
> https://www.animalgenome.org/cgi-bin/QTLdb/BT/download?d=0RxrNMbA8WzeDHyKIdCOL

This file need to be trimmed by removing innacurate regions (with start position higher than end and negative positions. We can use the following command:
>  awk '($2<$3&&$2>0){print $0}' QTL.bed > QTL_trimmed.bed

Necesitamos coger las DMCs unicamente con chr pos pos. Hay que filtrar el DMC delta FDR para que tenga esa pinta


Bedtools is used to find the overlapping regions using the intersect option:
> bedtools intersect -a QTL_trimmed.bed -b DMR.bed -wa | uniq > QTLsInDMR.tsv

We can also trim the file to only keep the coordinates and QTL associated to that region for further analyses:
> cut -f1,2,3,4 QTLsInDMR.tsv > QTLsInDMR_trimmed.tsv

## Enrichment

https://www.webgestalt.org/
