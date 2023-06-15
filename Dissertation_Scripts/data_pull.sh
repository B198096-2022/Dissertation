#This script will generate make a new directory for the incoming data files
#Download the specified data files
#Rename the data files 
#And generate a parameters file for the submit.jobs.sh script 
#The desired destination directory, files, and file names are 
#specified in the parameters file that is passed to submit_pull.sh


#!/bin/bash

#Establishing parameters 
while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

count=$((1))

#Loop through each of the SRR IDs as specified when running submit_pull.sh
while read id; do
  count2=${id: -2}
  #Make directory and set as working directory 
  mkdir /exports/eddie/scratch/s2249132/data/${big_dir}/${sub_dir}_${count}
  cd /exports/eddie/scratch/s2249132/data/${big_dir}/${sub_dir}_${count}
  #Download the files and rename them 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR_short}/0${count2}/${SRR_long}${id}/${SRR_long}${id}_1.fastq.gz
  mv ${SRR_long}${id}_1.fastq.gz ${SRR_long}${id}_S1_L001_R1_001.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR_short}/0${count2}/${SRR_long}${id}/${SRR_long}${id}_2.fastq.gz
  mv ${SRR_long}${id}_2.fastq.gz ${SRR_long}${id}_S1_L001_R2_001.fastq.gz
  #Move out one directory and make the parameters file for submit.jobs.sh
  cd /exports/eddie/scratch/s2249132/data/${big_dir}
  echo stage=0 > params.${sub_dir}_${count}
  echo tidy=0 >> params.${sub_dir}_${count}
  echo dir=/exports/eddie/scratch/s2249132/data/${big_dir} >> params.${sub_dir}_${count}
  echo project=${sub_dir}_${count} >> params.${sub_dir}_${count}
  echo n=0 >> params.${sub_dir}_${count}
  echo ref=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/references/refdata-gex-GRCh38-2020-A/ >> params.${sub_dir}_${count}
  echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> params.${sub_dir}_${count}
  echo destination=/exports/cmvm/datastore/scs/groups/pramacha-GROUP/Max_Hammer/cellranger/${big_dir} >> params.${sub_dir}_${count}
  echo srcdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/src >> params.${sub_dir}_${count}
  count=$((count+1))
done < ${sub_dir}_numbers.txt

