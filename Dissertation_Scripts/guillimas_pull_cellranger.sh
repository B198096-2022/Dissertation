#Guilliams

ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/049/SRR17375049/SRR17375049_1.fastq.gz


#This is the final command
#Run it in the scratch home space
#Mkdir the big_dir for all of the data first
#Both scripts used are in the data dir, but it is specified here
#The script, the project name, and the pull parameters file
./data/submit.pull.sh guilliams_healthy_cell pull_parameters.guilliams_healthy_cell

#Making the pull parameters file in home scratch dir
#big_dir is the dir for the full dataset
#sub_dir is for if there are multiple experimental groups (healthy vs fibrosis)
#SRR_short is the beginning of the SRR id that appears in the url
#SRR_long is the repeating portion of the full SRR id in the  url,
#       whatever is not specified by the project_numbers.txt
pull_parameters.guilliams_healthy_cell
echo big_dir=guilliams > pull_parameters.guilliams_healthy_cell
echo sub_dir=guilliams_healthy_cell >> pull_parameters.guilliams_healthy_cell
echo SRR_short=SRR173 >> pull_parameters.guilliams_healthy_cell
echo SRR_long=SRR173750 >> pull_parameters.guilliams_healthy_cell
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_healthy_cell

seq 32 33 >> guilliams_healthy_cell_numbers.txt
seq 46 47 >> guilliams_healthy_cell_numbers.txt


######################

./data/submit.pull.sh guilliams_fibrosis_cell pull_parameters.guilliams_fibrosis_cell

pull_parameters.guilliams_fibrosis_cell
echo big_dir=guilliams > pull_parameters.guilliams_fibrosis_cell
echo sub_dir=guilliams_fibrosis_cell >> pull_parameters.guilliams_fibrosis_cell
echo SRR_short=SRR173 >> pull_parameters.guilliams_fibrosis_cell
echo SRR_long=SRR173750 >> pull_parameters.guilliams_fibrosis_cell
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_fibrosis_cell

echo 38 > guilliams_fibrosis_cell_numbers.txt
seq 44 45 >> guilliams_fibrosis_cell_numbers.txt
seq 48 50 >> guilliams_fibrosis_cell_numbers.txt

######################

./data/submit.pull.sh guilliams_fibrosis_nuc pull_parameters.guilliams_fibrosis_nuc

pull_parameters.guilliams_fibrosis_nuc
echo big_dir=guilliams > pull_parameters.guilliams_fibrosis_nuc
echo sub_dir=guilliams_fibrosis_nuc >> pull_parameters.guilliams_fibrosis_nuc
echo SRR_short=SRR173 >> pull_parameters.guilliams_fibrosis_nuc
echo SRR_long=SRR173750 >> pull_parameters.guilliams_fibrosis_nuc
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_fibrosis_nuc

echo 11 > guilliams_fibrosis_nuc_numbers.txt
seq 54 58 >> guilliams_fibrosis_nuc_numbers.txt

######################

./data/submit.pull.sh guilliams_healthy_nuc pull_parameters.guilliams_healthy_nuc

pull_parameters.guilliams_fibrosis_nuc
echo big_dir=guilliams > pull_parameters.guilliams_healthy_nuc
echo sub_dir=guilliams_healthy_nuc >> pull_parameters.guilliams_healthy_nuc
echo SRR_short=SRR173 >> pull_parameters.guilliams_healthy_nuc
echo SRR_long=SRR173750 >> pull_parameters.guilliams_healthy_nuc
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_healthy_nuc

seq 51 53 >> guilliams_healthy_nuc_numbers.txt




######################

submit.pull.sh
#!/bin/bash
jid=$1
parameters=$2

while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

qsub -N pull.${jid} -cwd -l h_vmem=32g -V -v parameters=${parameters} -j y -o ${logdir}/ /exports/eddie/scratch/s2249132/data/data_pull.sh

######################
data_pull.sh
#!/bin/bash

while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

count=$((1))

while read id; do
  count2=${id: -2}
  mkdir /exports/eddie/scratch/s2249132/data/${big_dir}/${sub_dir}_${count}
  cd /exports/eddie/scratch/s2249132/data/${big_dir}/${sub_dir}_${count}
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR_short}/0${count2}/${SRR_long}${id}/${SRR_long}${id}_1.fastq.gz
  mv ${SRR_long}${id}_1.fastq.gz ${SRR_long}${id}_S1_L001_R1_001.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR_short}/0${count2}/${SRR_long}${id}/${SRR_long}${id}_2.fastq.gz
  mv ${SRR_long}${id}_2.fastq.gz ${SRR_long}${id}_S1_L001_R2_001.fastq.gz
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


#To run cell cellranger
#Move to the main scratch home
#Move param files and the project_files.txt to this space and then run
cp ./data/guilliams/params.* .
#
cp ./data/guilliams/params.guilliams_healthy_cell* .
seq 1 4 >> guilliams_healthy_cell_list.txt
while read numb; do echo guilliams_healthy_cell_${numb} >> guilliams_healthy_cell_files.txt; done < guilliams_healthy_cell_list.txt
while read sample; do ./src/submit.jobs.sh ${sample} params.${sample}; done < guilliams_healthy_cell_files.txt

#
cp ./data/guilliams/params.guilliams_fibrosis_cell* .
seq 1 6 >> guilliams_fibrosis_cell_list.txt
while read numb; do echo guilliams_fibrosis_cell_${numb} >> guilliams_fibrosis_cell_files.txt; done < guilliams_fibrosis_cell_list.txt
while read sample; do ./src/submit.jobs.sh ${sample} params.${sample}; done < guilliams_fibrosis_cell_files.txt

#
cp ./data/guilliams/params.guilliams_fibrosis_nuc* .
seq 1 6 >> guilliams_fibrosis_nuc_list.txt
while read numb; do echo guilliams_fibrosis_nuc_${numb} >> guilliams_fibrosis_nuc_files.txt; done < guilliams_fibrosis_nuc_list.txt
while read sample; do ./src/submit.jobs.sh ${sample} params.${sample}; done < guilliams_fibrosis_nuc_files.txt

#
cp ./data/guilliams/params.guilliams_healthy_nuc* .
seq 1 3 >> guilliams_healthy_nuc_list.txt
while read numb; do echo guilliams_healthy_nuc_${numb} >> guilliams_healthy_nuc_files.txt; done < guilliams_healthy_nuc_list.txt
while read sample; do ./src/submit.jobs.sh ${sample} params.${sample}; done < guilliams_healthy_nuc_files.txt
