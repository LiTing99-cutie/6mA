################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_Pacbio.Tetrahymena.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 18 Apr 2022 10:46:00 AM CST
################################################

#!/bin/sh 

set -eou pipefail

BASEMODS_PATH=/home/user/data2/lit/project/6mA/reproduce/smrt_output
SCRIPT=/home/user/data2/lit/project/6mA/reproduce/bin/binned_level_list_based.sh
BINNED_DIR=/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center

# ----

GENOME_NAME=T.thermophila
logDir=/home/user/data2/lit/project/6mA/log/iNPS/$GENOME_NAME
OUTPUT_DIR=/home/user/data/lit/project/6mA/feature/iNPS/$GENOME_NAME
[ -d $logDir ] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log

bash $SCRIPT \
-s "T.thermophila|" \
-o $OUTPUT_DIR \
-mA "$BASEMODS_PATH/SMRT-list-Tetrahymena.bed|" \
-fs /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
-fa /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta \
-ld $logDir \
-bd $BINNED_DIR \
> $logDir/run.log 2>&1

# nohup bash /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_list_based.T.thermophila.sh &