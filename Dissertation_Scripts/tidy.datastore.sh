#This script moves the outputs of cell ranger to the data store 
#And if specified will remove the remainig files from the scratch space 

#!/bin/bash

#Declare parameters passed as arguments 
while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

#If the destination doesn't exist, make a directory for it
if [ ! -e "${destination}" ]; then
  mkdir ${destination}
fi
#If the project directory doesn't exist, make it
if [ ! -e "${destination}/${project}" ]; then
  mkdir ${destination}/${project}
fi

# tidy is a variable in param file, do you want to delete in scratch after 
if [ ${tidy} -eq 1 ]; then
  rm -r ${dir}/${project}
fi

#Then this moves the the project directory from the scratch to the DataStore
rsync -avu ${dir}/${project}_output/* ${destination}/${project}
