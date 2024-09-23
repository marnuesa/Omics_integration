#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=30               # Number of cpus
#SBATCH -t 0-07:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=30G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o bed_create_%j.log # Standard output goes to this file

################################################################################
#
#               01-Bed_files_create.sh
#   This program is divide in three parts, that save files in order of priority:
#
#   The first part use the gff3 file of v4 Cucumis Melo genome, load in 
#   Cucurbit genomics, to generate bed files of differents features:
#
#   - Gene: Since start to end of the gene. But if the gene have 5'UTR or 
#   3'UTR it start at the end of 5'UTR and end at the start of 3'UTR
#   - Upstream: A 2000 bp region upstream of the gene start (if a 5' UTR is 
#   present, it is included in this region)
#   - Downstream: A 2000 bp region downstream of the gene end (if a 3' UTR 
#   is present, it is included in this region)
#
#   The second part, align all the Cucumir melo microRNA precursors in miRBase 
#   to the genome and generate bed files of three features of microRNA:
#
#   - microRNA: A region that starts where the precursor aligns and ends 
#   at the sum of the start position plus its length
#   - 5prime: A 1000 bp region upstream of the precursor start
#   - 3prime: A 1000 bp region downstream of the gene end 
#   
#   The third part, use a GTF of lncRNA in Cucumis melo to generate
#   a bed file with three features od lncRNA:
#
#   - lncRNA: A region that starts where the precursor aligns and ends 
#   at the sum of the start position plus its length
#   - 5primelncRNA: A 1000 bp region upstream of the precursor start
#   - 3primelncRNA: A 1000 bp region downstream of the gene end 
#   
#   At the end, the script add a column to all the files indicating the type of 
#   feature it is.
################################################################################


# GENE REGIONS
gene_file='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/tx2gene.txt'
gff_file='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/DHL92_v4.gff3'
bed_path='/home/nuezsal/Omics_integration/Integration_methylome/Results/Bed_files'
output_up="$bed_path/00-CMelon_DHL92_v4_upstream.bed"
output_down="$bed_path/06-CMelon_DHL92_v4_downstream.bed"
output_genes="$bed_path/01-CMelon_DHL92_v4_genes.bed"

mkdir -p $bed_path

# Clean output file if it exist
> "$output_up"
> "$output_down"
> "$output_genes"

# Read through the GFF file and process the lines
echo "...PROCESSING GFF3 FILE..."

cut -f 2 "$gene_file" | while read -r gene_id; do

	# Select only the rows of this gene
	grep "$gene_id" "$gff_file" > "temp_gene.gff3"
	lines=$(wc -l < "temp_gene.gff3")
	
	awk -F '\t' '{print NR, $1, $3, $4, $5, $7, $9}' "temp_gene.gff3" | while read line_num chr types start end strand gene_name; do
		# If gene is in forward strand
		if [ "$strand" == "+" ]; then
			
			# Control variable "three_prime_UTR"
			three_prime_UTR_processed=false

			if [ "$types" == "gene" ]; then
		    		
				# Coordinates of gene if there aren't 5' or 3' UTR regions
				start_gene=$start
				end_gene=$end
				
				# Extract up-stream regions 
		    		# Start of up-stream region
		    		if [ $start -gt 2000 ]; then
					start_upstream=$(($start - 2000))
				else
					start_upstream=0
				
				fi
				# End of upstream region
				end_upstream=$(($start - 1))
				
				# Extract down-stream regions 
		    		# Start of down-stream region
				start_downstream=$((end + 1))
			
				# End of down-stream region
				end_downstream=$(($end + 2000))
				
			elif [ "$types" == "five_prime_UTR" ]; then
				
				# Read the following lines to adjust the END of up-stream region if there is a five_prime_UTR or not
				end_upstream=$end
				
				# Change start gene if exist 5'UTR region
		        	if [ $end_upstream -ne $start ]; then
		        		start_gene=$(($end_upstream + 1))
		        	fi
					
			elif [ "$types" == "three_prime_UTR" ]; then
			
				if [ "$three_prime_UTR_processed" == false ]; then
			    		three_prime_UTR_processed=true
			    		
			    		# Read the following lines to adjust the start of down-stream region if there is a three_prime_UTR or not
			    		start_downstream=$start
			    		
					# Change end gene if exist 3'UTR region
		    			if [ $start_downstream -ne $end ]; then
		        			end_gene=$(($start_downstream - 1))
		        		fi
				fi
		    	fi
		        
		# If gene is in reverse strand
		else
			# Control variable "three_prime_UTR"
			five_prime_UTR_processed=false

			if [ "$types" == "gene" ]; then
		    		
				# Coordinates of gene if there aren't 5' or 3' UTR regions
				start_gene=$start
				end_gene=$end
				
				# Extract up-stream regions 
		    		# Start of up-stream region
		    		start_upstream=$(($end + 1))
				
				# End of upstream region
				end_upstream=$(($end + 2000))
				
				# Extract down-stream regions 
		    		# Start of down-stream region
		    		if [ $start -gt 2000 ]; then
					start_downstream=$(($start - 2000))
				else
					start_downstream=0
				fi
				
				# End of down-stream region
				end_downstream=$(($start - 1))
				
			elif [ "$types" == "five_prime_UTR" ]; then
			
				if [ "$five_prime_UTR_processed" == false ]; then
			    		five_prime_UTR_processed=true
			    		
			    		# Read the following lines to adjust the start of up-stream region if there is a five_prime_UTR or not
			    		start_upstream=$start
			    		
					# Change end gene if exist 5'UTR region
		    			if [ $start_upstream -ne $end ]; then
		        			end_gene=$(($start_upstream - 1))
		        		fi
				fi
					
			elif [ "$types" == "three_prime_UTR" ]; then
			
				# Read the following lines to adjust the END of down-stream region if there is a three_prime_UTR or not
				end_downstream=$end
				# Change start gene if exist 3'UTR region
		        	if [ $end_downstream -ne $start ]; then
		        		start_gene=$(($end_downstream + 1))
		        	fi				
		    	fi			
				
		fi
		
		if [ "$line_num" -eq "$lines" ]; then
			# Write gene at the files
			echo -e "$chr\t$start_upstream\t$end_upstream\t$gene_id" >> "$output_up"
			echo -e "$chr\t$start_downstream\t$end_downstream\t$gene_id" >> "$output_down"
			echo -e "$chr\t$start_gene\t$end_gene\t$gene_id" >> "$output_genes"
		
		fi
	done
	
done

echo "...END GENE FEATURE EXTRACTION..."

rm "temp_gene.gff3"

# MICRORNA REGIONS

Genome='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/'
hairpin_file='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/hairpin_cme_T.fa'
sam_file='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/microRNA_hairpin.sam'

output_5prime="$bed_path/02-microRNA_5prime.bed"
output_microRNA="$bed_path/03-microRNA_precursors.bed"
output_3prime="$bed_path/07-microRNA_3prime.bed"

# Alignment
# Bowtie Index
echo "...INDEXING C.MELON GENOME..."
bowtie-build $Genome/DHL92_genome_v4.fa $Genome/Melon_index

# Bowtie align
echo "...ALIGN MICRORNA PRECURSORS..."
bowtie -x $Genome/Melon_index --best -v 3 -k1 --no-unal -f $hairpin_file -S $sam_file

# Create bed file
# Clean output file if it exist
> "$output_5prime"
> "$output_microRNA"
> "$output_3prime"

# Read through the SAM file and process the line
echo "...EXTRACT MICRORNA FEATURES..."

awk -F '\t' '!/^@/ {print $1, $2, $3, $4, $10}' "$sam_file" | while read -r microRNA flag chr start seq; 
do
	if [ "$flag" == 0 ];then
		# Create variables to microRNAs
		start_forward=$start
		length=$(expr length "$seq")
		end_forward=$((start + $length))
		
		echo -e "$chr\t$start_forward\t$end_forward\t$microRNA" >> "$output_microRNA"
		
		# Create variables to 5' elements
		end_5=$(($start_forward - 1))
		if (($start_forward > 1000)); then
			start_5=$(($start_forward - 1000))
		else
			start_5=0
		fi
		
		echo -e "$chr\t$start_5\t$end_5\t$microRNA" >> "$output_5prime"
		
		# Create variables to 3' elements
		start_3=$(($end_forward + 1))
		end_3=$(($end_forward + 1000))
		
		echo -e "$chr\t$start_3\t$end_3\t$microRNA" >> "$output_3prime"
	
	elif [ "$flag" == 16 ]; then
	
		# Create variables to microRNAs
		length=$(expr length "$seq")
		start_reverse=$(($start - length))
		end_reverse=$start
		
		echo -e "$chr\t$start_reverse\t$end_reverse\t$microRNA" >> "$output_microRNA"
		
		# Create variables to 5' elements
		start_5=$(($end_reverse + 1000))
		end_5=$(($end_reverse + 1))
		
		echo -e "$chr\t$end_5\t$start_5\t$microRNA" >> "$output_5prime"
		
		# Create variables to 3' elements
		end_3=$(($start_reverse - 1))
		if (($start_reverse > 1000)); then
			start_3=$(($start_reverse - 1000))
		else
			start_3=0
		fi
		
		
		echo -e "$chr\t$start_3\t$end_3\t$microRNA" >> "$output_3prime"
	fi
done

echo "...END EXTRACT MICRORNA FEATURES..."

# LNCRNA REGIONS

gtf_file_lnc='/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/POTENTIAL_LNCRNAS_pred.gtf'
output_file="$bed_path/05-CMelon_DHL92_v4_lncRNA.bed"
prime5_file="$bed_path/04-CMelon_DHL92_v4_5primelncRNA.bed"
prime3_file="$bed_path/08-CMelon_DHL92_v4_3primelncRNA.bed"

# Clean output file if it exist
> "$output_file"
> "$prime5_file"
> "$prime3_file"

# Read through the GTF file and process the lines
echo "...EXTRACT LNCRNA FEATURES..."

awk -F '\t' '{print NR, $1, $3, $4, $5, $7, $9}' "$gtf_file_lnc" | while read line_num chr types start end strand lnc_name; do

	if [ "$strand" == "+" ]; then 	
		if [ "$types" == "transcript" ]; then
			
			# Extract the gene ID from the gene_name field
		    	lnc_id=$(echo "$lnc_name" | awk -F'[;]' '{print $1}' | awk -F' ' '{print $2}'| tr -d '"')	  		
			echo -e "$chr\t$start\t$end\t$lnc_id" >> "$output_file"
			
			# Extract upstream regions
			if (($start > 1000)); then
				start_up=$(($start - 1000))
			else
				start_up=0
			fi
			
			if (($start > 1)); then
				end_up=$(($start - 1))		
				echo -e "$chr\t$start_up\t$end_up\t$lnc_id" >> "$prime5_file"
			fi
			# Extract downstream regions
			start_down=$(($end + 1))
			end_down=$(($end + 1000))
			
			
			echo -e "$chr\t$start_down\t$end_down\t$lnc_id" >> "$prime3_file"		
		fi
			
	else
		if [ "$types" == "transcript" ]; then
			
			# Extract the gene ID from the gene_name field
		    	lnc_id=$(echo "$lnc_name" | awk -F'[;]' '{print $1}' | awk -F' ' '{print $2}'| tr -d '"')	  		
			echo -e "$chr\t$start\t$end\t$lnc_id" >> "$output_file"
			
			# Extract upstream regions
			start_up=$(($end + 1))
			end_up=$(($end + 1000))
			
			echo -e "$chr\t$start_up\t$end_up\t$lnc_id" >> "$prime5_file"
			
			# Extract downstream regions
			if (($start > 1000)); then
				start_down=$(($start - 1000))
			else
				start_down=0
			fi
			
			if (($start > 1)); then
				end_down=$(($start - 1))
				echo -e "$chr\t$start_down\t$end_down\t$lnc_id" >> "$prime3_file"	
			fi	
		fi
	
	fi		
done

echo "...END EXTRACT LNCRNA FEATURES..."

# MOVE RETROTRANSPOSON FILE

mv "/home/nuezsal/Omics_integration/Integration_methylome/Additional_info/CMelon_DHL92_v4_retrotransposons.bed" "$bed_path/09-CMelon_DHL92_v4_retrotransposons.bed"

# Add the column with the feature

echo "...ADD COLUMN WITH THE FEATURE..."

for file in "$bed_path"/*;
do
	echo "$file"	
	# Get the file name without the path
    	file_name=$(basename "$file")
    	
    	# Remove the .bed extension
	file_name_without_ext="${file_name%.bed}"

	# Split by _ and keep the last value
	name_out=$(echo "$file_name_without_ext" | rev | cut -d'_' -f1 | rev)
	echo "$name_out"
	awk -v OFS="\t" -v name="$name_out" '{print $0, name}' "$file" > "$bed_path/${file_name_without_ext}_new.bed"

	# Renombra el archivo temporal para reemplazar el original
	mv "$bed_path/${file_name_without_ext}_new.bed" "$file"
done


