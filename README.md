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
