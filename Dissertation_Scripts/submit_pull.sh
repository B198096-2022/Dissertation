#This script will submit the data pull commands to Eddie 

#!/bin/bash
jid=$1
parameters=$2

while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

qsub -N pull.${jid} -cwd -l h_vmem=32g -V -v parameters=${parameters} -j y -o ${logdir}/ /exports/eddie/scratch/s2249132/data/data_pull.sh
