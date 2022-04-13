################################################
#File Name: bin/run.motif.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 08 Apr 2022 08:15:42 PM CST
################################################

#!/bin/sh 

set -eou pipefail

script=/home/user/data2/lit/project/6mA/reproduce/bin/motif.sh

<<'!'
bash ${script} \
/home/user/data2/lit/project/6mA/reproduce/output/feature/iNPS/T.thermophila \
/home/user/data2/lit/project/6mA/reproduce/log/feature/iNPS/T.thermophila \
50 \
FALSE \
/home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta \
/home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
/home/user/data/lit/project/6mA/reproduce/iNPS_output \
T.thermophila

# nohup bash /home/user/data2/lit/project/6mA/reproduce/bin/run.motif.sh &
!

bash ${script} \
/home/user/data2/lit/project/6mA/reproduce/output/feature/iNPS/C.reinhardtii \
/home/user/data2/lit/project/6mA/reproduce/log/feature/iNPS/C.reinhardtii \
50 \
FALSE \
/home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa \
/home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
/home/user/data/lit/project/6mA/reproduce/iNPS_output \
C.reinhardtii

# nohup bash /home/user/data2/lit/project/6mA/reproduce/bin/run.motif.sh >log/run.motif.C.reinhardtii.log 2&>1 &
# time bash ${script}... >log...
# could be written better