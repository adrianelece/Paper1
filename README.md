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
