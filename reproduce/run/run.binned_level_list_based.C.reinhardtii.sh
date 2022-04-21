################################################
#File Name: /home/user/data2/lit/project/6mA/run/run.binned_level_Pacbio.Algae.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 15 Apr 2022 09:11:37 PM CST
################################################

#!/bin/sh 

set -eou pipefail

BASEMODS_PATH=/home/user/data2/lit/project/6mA/reproduce/smrt_output
SCRIPT=/home/user/data2/lit/project/6mA/reproduce/bin/binned_level_list_based.sh
BINNED_DIR=/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center


# ----

GENOME_NAME=C.reinhardtii
logDir=/home/user/data2/lit/project/6mA/log/iNPS/$GENOME_NAME
OUTPUT_DIR=/home/user/data/lit/project/6mA/feature/iNPS/$GENOME_NAME
[ -d $logDir ] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log


bash $SCRIPT \
-s "C.reinhardtii|" \
-o $OUTPUT_DIR \
-mA "$BASEMODS_PATH/SMRT-list-Algae.bed|" \
-fs /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
-fa /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa \
-ld $logDir \
-bd $BINNED_DIR \
> $logDir/run.log 2>&1

# nohup bash /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_list_based.C.reinhardtii.sh &