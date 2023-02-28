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
