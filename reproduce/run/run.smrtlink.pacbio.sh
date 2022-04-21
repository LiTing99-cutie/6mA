################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/bin/run.smrtlink.pacbio.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 06 Apr 2022 04:16:50 PM CST
################################################

#!/bin/sh 


set -eou pipefail

script=/home/user/data2/lit/project/6mA/reproduce/bin/smrtlink.pacbio.sh

bash ${script} \
/home/user/data/lit/project/6mA/reproduce/rawdata/Pacbio/bam/Algae.subreads.bam \
/home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa \
/home/user/data2/lit/project/6mA/reproduce/output/Pacbio/Algae/Algae.ccs.fasta \
/home/user/data2/lit/project/6mA/reproduce/output/Pacbio/ \
/home/user/data2/lit/project/6mA/reproduce/log

# nohup bash /home/user/data2/lit/project/6mA/reproduce/bin/run.smrtlink.pacbio.sh &