#!/bin/sh 
################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_Pacbio.SMRT-1.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 16 Apr 2022 02:54:31 PM CST
################################################

<<'!'
set -eou pipefail

BASEMODS_PATH=/home/user/data2/lit/project/6mA/reproduce/smrt_output
SCRIPT=/home/user/data2/lit/project/6mA/reproduce/bin/binned_level_list_based.sh
BINNED_DIR=/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center

# ----

GENOME_NAME=C.elegans
mA_list_PATH=/home/user/data2/lit/project/6mA/mA_list
logDir=/home/user/data2/lit/project/6mA/log/danpos/normed/$GENOME_NAME
OUTPUT_DIR=/home/user/data/lit/project/6mA/feature/danpos/normed/$GENOME_NAME

[ -d $logDir ] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log

# -s Mnase sample
bash $SCRIPT \
-s $(for i in `seq 1 2`;do printf "PC10_LY_MCC-M$i|";done) \
-o $OUTPUT_DIR \
-mA "$BASEMODS_PATH/SMRT-list-SMRT-1.bed|$mA_list_PATH/sig_gt2_fc0_wt1-wt2_merge_correct-strand-1.bed.input|$mA_list_PATH/sig_gt2_fcnot0_wt1-wt2_merge_correct-strand-2.bed.input|$mA_list_PATH/merge_correct-strand-3.bed.input" \
-fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
-fa /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa \
-ld $logDir \
-bd $BINNED_DIR \
> $logDir/run.log 2>&1

# nohup bash /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_list_based.C.elegans.sh &
!

set -eou pipefail

BASEMODS_PATH=/home/user/data2/lit/project/6mA/reproduce/smrt_output
SCRIPT=/home/user/data2/lit/project/6mA/reproduce/bin/binned_level_list_based.sh
GENOME_NAME=C.elegans
mA_list_PATH=/home/user/data2/lit/project/6mA/mA_list

for bin_number in 10 8 6 4;do
    OUTPUT_DIR=/home/user/data/lit/project/6mA/feature/danpos/normed/${bin_number}bin/$GENOME_NAME
    logDir=/home/user/data2/lit/project/6mA/log/danpos/normed/${bin_number}bin/$GENOME_NAME
    BINNED_DIR=/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin
    [ -d $logDir ] || mkdir -p $logDir
    [[ -f $logDir/run.log ]] && rm -rf $logDir/run.log
    time 
    ( bash $SCRIPT \
    -s $(for i in $(seq 1 2);do printf "PC10_LY_MCC-M$i|";done) \
    -o $OUTPUT_DIR \
    -mA "$BASEMODS_PATH/SMRT-list-SMRT-1.bed|$mA_list_PATH/sig_gt2_fc0_wt1-wt2_merge_correct-strand-1.bed.input|$mA_list_PATH/sig_gt2_fcnot0_wt1-wt2_merge_correct-strand-2.bed.input|$mA_list_PATH/sig_gt2_fcnot0_wt1-wt2_merge_correct-strand-2-shuf.bed.input.sorted|$mA_list_PATH/merge_correct-strand-3.bed.input" \
    -fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
    -fa /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa \
    -ld $logDir \
    -bd $BINNED_DIR \
    -n $bin_number \
    > $logDir/run.log 2>&1 ) &
done

wait

# nohup bash /home/user/data2/lit/project/6mA/reproduce/run/run.binned_level_list_based.C.elegans.sh &
