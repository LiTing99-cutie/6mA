################################################
#File Name: bin/annotation_download.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 07 Mar 2022 07:43:35 PM CST
################################################

#!/bin/sh 

echo `date`

processed_data_dir=/home/user/data/lit/project/6mA/reproduce/processed_data/

gtf_dir=/home/user/data/lit/database/public/annotation/gtf

genome_data_dir_c=/home/user/data/lit/database/public/genome/C.reinhardtii
genome_data_dir_t=/home/user/data/lit/database/public/genome/T.thermophila
### download annotation

# T. thermophila
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014_assembly.fasta -P $genome_data_dir_t &
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014.gff3 -P $gtf_dir &

# C. reinhardtii

# download in windows

### download code and processed data

wget \
https://zenodo.org/record/5838427/files/Code_and_processed_data_for_Critical_assessment_of_DNA_adenine_methylation_in_eukaryotes_using_quantitatived_econvolution_new_license.zip?download=1https://zenodo.org/record/5838427/files/Code_and_processed_data_for_Critical_assessment_of_DNA_adenine_methylation_in_eukaryotes_using_quantitatived_econvolution_new_license.zip?download=1 \
-O $processed_data_dir/processed_data.zip &

wait

echo "done"

# nohup bash bin/annotation_download.sh > log/annotation_download.log 2>&1 &