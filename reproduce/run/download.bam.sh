#!/bin/sh 

output_dir=/home/user/data2/lit/DATA/project/6mA/reproduce/rawdata/Pacbio/

# test
# time (prefetch SRR15679411 -O /home/user/data2/lit/DATA/project/6mA/reproduce/rawdata/Pacbio/) > prefetch.log 2>&1 & 

# time (prefetch SRR15687574 --max-size 40G -O $output_dir)

time ( for SRR in SRR15679411 SRR15687574;do
    sam-dump $output_dir/$SRR/$SRR.sra > $output_dir/$SRR/$SRR.sam
    # samtools view -buSH $output_dir/$SRR/$SRR.sam > $output_dir/$SRR/$SRR.bam
done )

# nohup bash bin/download.bam.sh >> log/download.bam.log 2>&1 &

# line 10-13
# nohup bash bin/download.bam.sh >> log/download.bam.log 2>&1 & 

mkdir /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/SRR15687574
mkdir /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/SRR15679411 
nohup aws s3 sync s3://datanano/SRR15687574/ /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/SRR15687574 2> log/aws_sync.SRR15687574.err.log 1> log/aws_sync.SRR15687574.log &
nohup aws s3 sync s3://datanano/SRR15679411/ /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/SRR15679411 2> log/aws_sync.SRR15679411.err.log 1> log/aws_sync.SRR15679411.log &
