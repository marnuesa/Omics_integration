#!/bin/bash

####### FUNCTIONS

task_subset(){
	echo -e "\t"$1"..."
	zcat $2/02-Clean_data/$1\_1.trimmed.fq.gz | head -n 8000000 | gzip > $2/SUBSET/$1\_1.trimmed.fq.gz
	zcat $2/02-Clean_data/$1\_2.trimmed.fq.gz | head -n 8000000 | gzip > $2/SUBSET/$1\_2.trimmed.fq.gz
}

"$@"

