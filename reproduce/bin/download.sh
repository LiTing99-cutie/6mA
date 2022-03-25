################################################
#File Name: bin/download.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 07 Mar 2022 05:03:26 PM CST
################################################

#!/bin/sh 

HOME=/home/user/BGM/lixs
rawdata_dir=/home/user/data/lit/project/6mA/reproduce/rawdata/Mnase

<<'!'
pushd $rawdata_dir

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR199/009/SRR1994849/SRR1994849_1.fastq.gz . && mv SRR1994849_1.fastq.gz SRR1994849_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR199/009/SRR1994849/SRR1994849_2.fastq.gz . && mv SRR1994849_2.fastq.gz SRR1994849_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR199/000/SRR1994850/SRR1994850_1.fastq.gz . && mv SRR1994850_1.fastq.gz SRR1994850_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR199/000/SRR1994850/SRR1994850_2.fastq.gz . && mv SRR1994850_2.fastq.gz SRR1994850_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR533/002/SRR5337752/SRR5337752_1.fastq.gz . && mv SRR5337752_1.fastq.gz SRR5337752_GSM2534785_MNase-seq_Tetrahymena_thermophila_MNase-Seq_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR533/002/SRR5337752/SRR5337752_2.fastq.gz . && mv SRR5337752_2.fastq.gz SRR5337752_GSM2534785_MNase-seq_Tetrahymena_thermophila_MNase-Seq_2.fastq.gz

popd
!

# nohup bash bin/download.sh >> log/download.log 2>&1 &

pushd $rawdata_dir

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR156/011/SRR15679411/SRR15679411_subreads.fastq.gz . && mv SRR15679411_subreads.fastq.gz SRR15679411_Pacbio_SEQ2_Short_insert_library_Chlamydomonas_reinhardtii_CC-125_mt__subreads.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR156/074/SRR15687574/SRR15687574_subreads.fastq.gz . && mv SRR15687574_subreads.fastq.gz SRR15687574_Pacbio_SEQ2_Short_insert_library_Tetrahymena_thermophila_CU428_subreads.fastq.gz

popd

# nohup bash bin/download.sh >> log/download.log 2>&1 &

mkdir /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio

nohup cp /home/user/data/lit/project/6mA/reproduce/rawdata/Mnase/*Pacbio* /home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio &

wait

rm -rf /home/user/data/lit/project/6mA/reproduce/rawdata/Mnase/*Pacbio*