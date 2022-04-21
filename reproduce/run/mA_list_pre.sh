################################################
#File Name: reproduce/run/mA_list_pre.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 19 Apr 2022 02:27:03 PM CST
################################################

#!/bin/sh 

set -eou pipefail


# SMRT

for sample in SMRT-1 Algae Tetrahymena;do
# 6mA sample
BASEMODS_PATH=/home/user/data2/lit/project/6mA/reproduce/smrt_output
BASEMODS=$BASEMODS_PATH/$sample.basemods.gff
less $BASEMODS | grep -v '#' | \
awk -F ';' -v OFS='\t' '{print $1,$3}' | sed 's/IPDRatio=//g' | \
awk -v OFS='\t' '{if ($3=="m6A") print $1,$4-1,$4,$10}' | sort -k1,1 -k2,2n > $BASEMODS_PATH/SMRT-list-$sample.bed
done 

# in_house final mA list

for mA_list in `ls /home/user/data2/lit/project/6mA/mA_list/*strand*bed`;do
less $mA_list | sed 's/chr//g;s/M/MtDNA/g' | awk -v OFS='\t' '{print $1,$2,$3,($5+$8)/($4+$7)}' | \
sort -k1,1 -k2,2n > $mA_list.input
done
