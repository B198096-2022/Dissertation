#This script just copies files from the data store to the scratch space
#It was not used for this project but is available as a functionality of submit.jobs.sh

#!/bin/bash

while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

if [ ! -e "${dir}/${project}" ]; then
  mkdir ${dir}/${project}
fi

#Is source supposed to be dir ...? It is where the fastqc files are kept 
rsync -vu ${source}/*.fastq.gz ${dir}/${project}
