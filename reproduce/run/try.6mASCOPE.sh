

#!/bin/sh

# wd：/home/user/data2/lit/

singularity shell --bind $PWD /home/user/data2/lit/6mASCOPE/

module load singularity/3.6.4    # Required only singularity/3.6.4 is a dynamic environment module. 
singularity pull 6mASCOPE.sif library://fanglabcode/default/6mascope:latest    # Download the image from cloud.sylabs.io; Make sure you have the network connection
singularity build --sandbox 6mASCOPE 6mASCOPE.sif     # Create a writable container named 6mASCOPE
singularity run --no-home -w /home/user/data2/lit/6mASCOPE    # Start an interactive shell to use 6mASCOPE, type `exit` to leave
# download error
init_6mASCOPE    #Inside the container; Only required once when start using 6mASCOPE
source run_6mASCOPE    #Inside the container; Required every time when running 6mASCOPE

singularity build --sandbox 6mASCOPE_s 6mASCOPE.sif 
singularity run --no-home -w 6mASCOPE_s

source run_6mASCOPE 
# Usage: 6mASCOPE <subtasks> [parameters]
#     subtasks:
#       - deplex: De-multiplex subreads into individal files.
#       - ccs: Generate CCS .bam and .fasta files from short insert library sequencing data.
#       - ipd: Summarize the IPD and QV information for each molecule by reference-free method.
#       - contam: Estimate the contamination in the current sample.
#       - quant: Quantify deconvolution of the sample with defined groups.

6mASCOPE get_test_data

# GnuTLS: The TLS connection was non-properly terminated.
# Unable to establish SSL connection.

6mASCOPE contam -c test.ccs.fasta -r test.ref.fasta -o test.contam.estimate.txt

6mASCOPE quant -c test.ccs.fasta -i test.IPD.out.A -o test -r test.ref.fasta -s subgroup.txt

6mASCOPE contam -c ./3_test-data/test.ccs.fasta.gz -r ./3_test-data/test.ref.fasta.gz -o test.contam.estimate.txt
6mASCOPE contam -c /home/user/data2/lit/DATA/project/6mA/reproduce/processed_data/submitted_scripts_github2/3_test-data/test.ccs.fasta \
-r /home/user/data2/lit/DATA/project/6mA/reproduce/processed_data/submitted_scripts_github2/3_test-data/test.ref.fasta \
-o test.contam.estimate.txt

gunzip -c  3_test-data/test.ccs.fasta.gz >3_test-data/test.ccs.fasta
gunzip -c  3_test-data/test.ref.fasta.gz >3_test-data/test.ref.fasta


6mASCOPE contam -c /home/6mASCOPE/test/test.ccs.fasta \
-r /home/6mASCOPE/test/test.ref.fasta \
-o test.contam.estimate.txt


dir_1=/home/user/data/lit/project/6mA/reproduce/processed_data/submitted_scripts_github2/1_fanglab_6mASCOPE_github/SLpackage
dir_2=/home/user/data/lit/project/6mA/reproduce/processed_data/submitted_scripts_github2/4_nt-pre-build-tar

dir_3=/home/user/data2/lit/6mASCOPE/home/6mASCOPE/SLpackage
dir_4=/home/user/data2/lit/6mASCOPE/home/6mASCOPE/database/nt-pre-build

nohup cp -r $dir_1/* $dir_3 &
nohup cp -r $dir_2/* $dir_4 &

for dir in $dir_1 $dir_2 $dir_3 $dir_4;do
nohup md5sum.sh $dir &
done


singularity shell --bind $PWD,/home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/:/bam /home/user/data2/lit/6mASCOPE/ 
source run_6mASCOPE 

nohup  md5sum.sh $dir_2 $dir_4 &

# run error
singularity run --no-home -w --bind /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/:/bam /home/user/data2/lit/6mASCOPE/ 

# not work
singularity run --no-home /home/user/data2/lit/6mASCOPE/ 

# not work
singularity run -w /home/user/data2/lit/6mASCOPE/ 

# work
singularity run --no-home -w /home/user/data2/lit/6mASCOPE/ 

./6mASCOPE_s/.singularity.d

./tmp_1/6mASCOPE/.singularity.d

#! /bin/bash
source /home/6mASCOPE/.bashrc 
eval "$(conda shell.bash hook)"
conda activate 6mASCOPE
6mASCOPE $@

nohup code/ipd.sh -r /home/user/data2/lit/DATA/project/6mA/reproduce/rawdata/Pacbio/bam/Algae.subreads.bam -f /home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/Algae.ccs.fasta -o Algae.Ipd.out -p 40 >ipd.Algae.log 2>&1 & 