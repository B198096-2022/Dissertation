#!/bin/bash

#This is the first job that coordinates everyting else
#You submit <jobid> and <parameters.file>

#The parameters.file has a list of varaibles that will be used for the rest of the pipeline

#-lt means less than, so if you have less than 2 arguments it says you are missing them
if [ $# -lt 2 ]; then
  echo "Invalid arguments supplied"
  echo "Script usage: 'submit.jobs.sh <job.id> <parameters.file>'"
  exit 0
else
  #Defining the job id is just the name that Eddie disaplays
  #parameters contains the variables
  jid=$1
  parameters=$2
fi

#-a specfies an array
#p is temporarily the word "stage" and next line assigns stage to being 0 (value) and stage = 0 is kept
#And we rest p, do the next one and go on maintaining variable names, p is just temporary
while IFS='=' read -a pArray; do
  p=${pArray[0]}
  declare ${p}="${pArray[1]}"
done < ${parameters}

if [ ${stage} -eq 1 ]; then
  #-N job name -cwd -q is defining which que (staging, so can see DS and eddie) -V inherits the local environment
  #-v defines variable names which are parameters
  #-j is "join" to join error and output files so it is all in one file
  #stage variable we set, is the file already in staging
  #srcdir is where the scripts are
  #So this if the files aren't already in staging, move them
  qsub -N cp.${jid} -cwd -q staging -V -v parameters=${parameters} -j y -o ${logdir}/ ${srcdir}/cp.datastore.sh
fi
#-gold_jid is looking at past job and waiting until it is done to run
#-l is resources allocation, so it is the limit of the job (ensure high memory node so you don't max it out)
qsub -N call.${jid} -hold_jid cp.${jid} -cwd -l h_vmem=32g -V -v parameters=${parameters} -j y -o ${logdir}/ ${srcdir}/call.sh
qsub -N tidy.${jid} -hold_jid call.${jid} -cwd -q staging -V -v parameters=${parameters} -j y -o ${logdir}/ ${srcdir}/tidy.datastore.sh
