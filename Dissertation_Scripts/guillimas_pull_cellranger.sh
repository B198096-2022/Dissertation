#Guilliams


############################################################
#                                                          #
#        Now run the submit.jobs.sh script                 #
#                                                          #
############################################################

#This is the final command
#Run it in the scratch home space
#Mkdir the big_dir for all of the data first
#Both scripts used are in the data dir, but it is specified here
#The script, the project name, and the pull parameters file

#./data/submit.pull.sh guilliams_healthy_cell pull_parameters.guilliams_healthy_cell

#Making the pull parameters file in home scratch dir
#big_dir is the dir for the full dataset
#sub_dir is for if there are multiple experimental groups (healthy vs fibrosis)
#SRR_short is the beginning of the SRR id that appears in the url
#SRR_long is the repeating portion of the full SRR id in the  url,
#       whatever is not specified by the project_numbers.txt

##########################
#  Healthy Cell samples  #
##########################
#pull_parameters.guilliams_healthy_cell

echo big_dir=guilliams > pull_parameters.guilliams_healthy_cell
echo sub_dir=guilliams_healthy_cell >> pull_parameters.guilliams_healthy_cell
echo SRR_short=SRR173 >> pull_parameters.guilliams_healthy_cell
echo SRR_long=SRR173750 >> pull_parameters.guilliams_healthy_cell
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_healthy_cell

#Specify the ends of the SRR accessions on ENA 

seq 32 33 >> guilliams_healthy_cell_numbers.txt
seq 46 47 >> guilliams_healthy_cell_numbers.txt

./data/submit.pull.sh guilliams_healthy_cell pull_parameters.guilliams_healthy_cell

######################

###########################
#  Fibrosis cell samples  #
###########################


#pull_parameters.guilliams_fibrosis_cell
echo big_dir=guilliams > pull_parameters.guilliams_fibrosis_cell
echo sub_dir=guilliams_fibrosis_cell >> pull_parameters.guilliams_fibrosis_cell
echo SRR_short=SRR173 >> pull_parameters.guilliams_fibrosis_cell
echo SRR_long=SRR173750 >> pull_parameters.guilliams_fibrosis_cell
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_fibrosis_cell

echo 38 > guilliams_fibrosis_cell_numbers.txt
seq 44 45 >> guilliams_fibrosis_cell_numbers.txt
seq 48 50 >> guilliams_fibrosis_cell_numbers.txt

./data/submit.pull.sh guilliams_fibrosis_cell pull_parameters.guilliams_fibrosis_cell

######################

###########################
#  Healthy Nuc  samples   #
###########################


#pull_parameters.guilliams_fibrosis_nuc
echo big_dir=guilliams > pull_parameters.guilliams_fibrosis_nuc
echo sub_dir=guilliams_fibrosis_nuc >> pull_parameters.guilliams_fibrosis_nuc
echo SRR_short=SRR173 >> pull_parameters.guilliams_fibrosis_nuc
echo SRR_long=SRR173750 >> pull_parameters.guilliams_fibrosis_nuc
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_fibrosis_nuc

echo 11 > guilliams_fibrosis_nuc_numbers.txt
seq 54 58 >> guilliams_fibrosis_nuc_numbers.txt

./data/submit.pull.sh guilliams_fibrosis_nuc pull_parameters.guilliams_fibrosis_nuc

######################

###########################
#  Fibrosis Nuc  samples  #
###########################


#pull_parameters.guilliams_fibrosis_nuc
echo big_dir=guilliams > pull_parameters.guilliams_healthy_nuc
echo sub_dir=guilliams_healthy_nuc >> pull_parameters.guilliams_healthy_nuc
echo SRR_short=SRR173 >> pull_parameters.guilliams_healthy_nuc
echo SRR_long=SRR173750 >> pull_parameters.guilliams_healthy_nuc
echo logdir=/exports/cmvm/eddie/scs/groups/pramacha-GROUP/max/cellranger/logs >> pull_parameters.guilliams_healthy_nuc

seq 51 53 >> guilliams_healthy_nuc_numbers.txt

./data/submit.pull.sh guilliams_healthy_nuc pull_parameters.guilliams_healthy_nuc

######################################################################
#                                                                    #
#              Now run the submit.jobs.sh script                     #
#                                                                    #
# Requires specifying the .src/ directory where scripts are located  #
#                                                                    #
######################################################################

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
