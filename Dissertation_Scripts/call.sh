#This script calls the cellranger function 
#to align fastq files to a reference genome
#And outputs count matrices 

#!/bin/bash

#Checks if the passed arguments contain all of the variables it needs
#pull in the arguments and declare them 
while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

#Load in cellranger
module load igmm/apps/cellranger/7.0.0

#specify sample labels if there are mutliple samples in teh fastq directory 
if [ ${n} -gt 0 ]; then
  samples="${pref}-1"; for i in `seq 2 ${n}`; do samples="${samples},${pref}-${n}"; done
#Then run cell ranger count
  cd ${dir}
  cellranger count --id=${project}_output \
    --transcriptome=${ref} \
    --fastqs=${dir}/${project} \
    --sample=${samples}
#Run cell ranger if there is a single sample in the fastq directory 
else
  cd ${dir}
  cellranger count --id=${project}_output \
    --transcriptome=${ref} \
    --fastqs=${dir}/${project}
fi
