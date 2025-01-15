#!/bin/bash

################################################################################
#
#               01-Bed_files_create.sh
#
#   Este script genera archivos BED a partir de un archivo GFF3, permitiendo 
#   especificar las rutas de entrada y salida mediante argumentos.
#
################################################################################

# Función de ayuda
usage() {
  echo "Uso: $0 -g <gene_file> -f <gff_file> -o <output_directory>"
  echo
  echo "  -g  Ruta al archivo tx2gene.txt"
  echo "  -f  Ruta al archivo gff"
  echo "  -o  Directorio de salida para los archivos BED"
  echo "  -h  Mostrar esta ayuda"
  exit 1
}

# Procesar argumentos
while getopts ":g:f:o:h" opt; do
  case ${opt} in
    g ) gene_file=$OPTARG ;;
    f ) gff_file=$OPTARG ;;
    o ) bed_path=$OPTARG ;;
    h ) usage ;;
    \? ) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    : ) echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done

# Verificar que los argumentos requeridos fueron proporcionados
if [[ -z "$gene_file" || -z "$gff_file" || -z "$bed_path" ]]; then
  echo "Error: Faltan argumentos obligatorios." >&2
  usage
fi

# Definir rutas de salida
output_up="$bed_path/CMelon_DHL92_v4_upstream.bed"
output_genes="$bed_path/CMelon_DHL92_v4_genes.bed"

mkdir -p "$bed_path"

# Limpiar archivos de salida si existen
> "$output_up"
> "$output_genes"

# Procesar archivo GFF3
echo "...Procesando archivo GFF3..."

cut -f 2 "$gene_file" | while read -r gene_id; do
  grep "$gene_id" "$gff_file" > "temp_gene.gff3"
  lines=$(wc -l < "temp_gene.gff3")

  awk -F '\t' '{print NR, $1, $3, $4, $5, $7, $9}' "temp_gene.gff3" | while read line_num chr types start end strand gene_name; do
    if [ "$strand" == "+" ]; then
      three_prime_UTR_processed=false

      if [ "$types" == "gene" ]; then
        start_gene=$start
        end_gene=$end
        start_upstream=$((start > 2000 ? start - 2000 : 1))
        end_upstream=$((start - 1))
      elif [ "$types" == "five_prime_UTR" ]; then
        end_upstream=$end
        start_gene=$((end_upstream + 1))
      elif [ "$types" == "three_prime_UTR" ] && [ "$three_prime_UTR_processed" == false ]; then
        three_prime_UTR_processed=true
        start_downstream=$start
        end_gene=$((start_downstream - 1))
      fi
    else
      five_prime_UTR_processed=false

      if [ "$types" == "gene" ]; then
        start_gene=$start
        end_gene=$end
        start_upstream=$((end + 1))
        end_upstream=$((end + 2000))
      elif [ "$types" == "five_prime_UTR" ] && [ "$five_prime_UTR_processed" == false ]; then
        five_prime_UTR_processed=true
        start_upstream=$start
        end_gene=$((start_upstream - 1))
      elif [ "$types" == "three_prime_UTR" ]; then
        end_downstream=$end
        start_gene=$((end_downstream + 1))
      fi
    fi

    if [ "$line_num" -eq "$lines" ]; then
      echo -e "$chr\t$((start_upstream - 1))\t$((end_upstream - 1))\t$strand\t$gene_id" >> "$output_up"
      echo -e "$chr\t$((start_gene - 1))\t$((end_gene - 1))\t$strand\t$gene_id" >> "$output_genes"
    fi
  done
done

echo "...Extracción completada..."

rm "temp_gene.gff3"
