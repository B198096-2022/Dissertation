#!/bin/bash



#Does this contain all of the variables we want? samples, project name, destination, etc?
#pull in the arguments?
while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

#Load in cellranger
module load igmm/apps/cellranger/7.0.0

#Is this then making sample labels?
if [ ${n} -gt 0 ]; then
  samples="${pref}-1"; for i in `seq 2 ${n}`; do samples="${samples},${pref}-${n}"; done
#Then run cell ranger count
  cd ${dir}
  cellranger count --id=${project}_output \
    --transcriptome=${ref} \
    --fastqs=${dir}/${project} \
    --sample=${samples}

else
  cd ${dir}
  cellranger count --id=${project}_output \
    --transcriptome=${ref} \
    --fastqs=${dir}/${project}
fi
