#!/bin/bash

####### FUNCTIONS

task_Genome_preparation(){
	## Remove the previous stdout log file.
	if [ -f "$2/Logs/Genome_preparation_stdout.log" ]; then
	    rm $2/Logs/Genome_preparation_stdout.log
	fi
	## Copy the genome file to the Genome preparation folder.
	cp $1/Genome/cme.fa $2/01-Genome_preparation/
	## Prepare genome. Parallel 1 uses really 2 cores, so parallel 8 uses 16 cores.
	bismark_genome_preparation --parallel 8 --verbose $2/01-Genome_preparation/ >> $2/Logs/Genome_preparation_stdout.log 2>&1
}

task_Alignment_paired_end(){
	echo -e "\t"$1"..."
	## Create the temporal folder.
	if [ -d "$3/02-Alignment/$1" ]; then
		rm -r $3/02-Alignment/$1
	fi
	mkdir $3/02-Alignment/$1
	## Remove the previous stdout log file.
	if [ -f "$3/Logs/Alignment_$1_stdout.log" ]; then
	    rm $3/Logs/Alignment_$1\_stdout.log
	fi
	## Bismark. Parallel 1 uses ~5 cores and ~10Gb, so parallel 4 uses ~20 cores and ~40Gb.
	bismark \
		--genome_folder $3/01-Genome_preparation/ \
		-1 $2/SUBSET/$1\_1.trimmed.fq.gz \
		-2 $2/SUBSET/$1\_2.trimmed.fq.gz \
		--local \
		--output_dir $3/02-Alignment/ \
		--temp_dir $3/02-Alignment/$1 \
		--parallel 4 >> $3/Logs/Alignment_$1\_stdout.log 2>&1
	## Rename.
	mv $3/02-Alignment/$1\_1.trimmed_bismark_bt2_pe.bam $3/02-Alignment/$1.trimmed.bismark_bt2_pe.bam
	mv $3/02-Alignment/$1\_1.trimmed_bismark_bt2_PE_report.txt $3/02-Alignment/$1.trimmed.bismark_bt2_PE_report.txt
	## Remove temporal folder.
	rm -r $3/02-Alignment/$1
}

task_Deduplication_paired_end(){
	echo -e "\t"$1"..."
	## Remove the previous stdout log file.
	if [ -f "$2/Logs/Deduplication_$1_stdout.log" ]; then
	    rm $2/Logs/Deduplication_$1\_stdout.log
	fi
	## Deduplication.
	deduplicate_bismark --paired --bam $2/02-Alignment/$1.trimmed.bismark_bt2_pe.bam --output_dir $2/03-Deduplication/ >> $2/Logs/Deduplication_$1\_stdout.log 2>&1
}

task_Methylation_extractor_paired_end(){
	echo -e "\t"$1"..."
	## Create the temporal folder.
	if [ -d "$2/04-Methylation_extractor/$1" ]; then
		rm -r $2/04-Methylation_extractor/$1
	fi
	mkdir $2/04-Methylation_extractor/$1
	## Remove the previous stdout log file.
	if [ -f "$2/Logs/Methylation_extractor_$1_stdout.log" ]; then
	    rm $2/Logs/Methylation_extractor_$1\_stdout.log
	fi
	## Methylation extraction. Parallel 1 uses really 3 cores, so parallel 8 uses 24 cores.
	bismark_methylation_extractor --parallel 8 --paired-end --comprehensive --no_overlap --ignore_3prime_r2 2 --cytosine_report --bedGraph --CX $2/03-Deduplication/$1.trimmed.bismark_bt2_pe.deduplicated.bam --genome_folder $2/01-Genome_preparation/ --output $2/04-Methylation_extractor/$1/ >> $2/Logs/Methylation_extractor_$1\_stdout.log 2>&1
	cp $2/04-Methylation_extractor/$1/* $2/04-Methylation_extractor/
	## Remove temporal folder.
	rm -r $2/04-Methylation_extractor/$1
}

task_QC(){
	echo -e "\t"$1"..."
	## Create the temporal folder.
	if [ -d "$2/05-Global_reports/$1" ]; then
		rm -r $2/05-Global_reports/$1
	fi
	mkdir $2/05-Global_reports/$1
	## Remove the previous stdout log file.
	if [ -f "$2/Logs/QC_$1_stdout.log" ]; then
	    rm $2/Logs/QC_$1\_stdout.log
	fi
	## Individual QC (by sample).
	cd $2/05-Global_reports/$1/
	cp $2/02-Alignment/*E_report.txt ./
	cp $2/03-Deduplication/**.deduplication_report.txt ./
	cp $2/04-Methylation_extractor/*M-bias.txt ./
	cp $2/04-Methylation_extractor/*_splitting_report.txt ./
	bismark2report >> $2/Logs/QC_$1\_stdout.log 2>&1
	cd ..
	cp ./$1/* ./
	rm -r $1
}

"$@"

