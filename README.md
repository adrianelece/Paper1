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
python ../DatabyDir.py -r /datasets/work/af-dairy-cpg/work/INIA/reanalysesMarch/7x/allregions.txt -s /datasets/work/af-dairy-cpg/work/INIA/reanalysesMarch/7x/ -p "network*.tsv" -o /datasets/work/af-dairy-cpg/work/INIA/reanalysesMarch/7x/AllSamples.tsv
```
