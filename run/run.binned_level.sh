################################################
#File Name: /home/user/data2/lit/project/6mA/bin/run.binned_level.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 12 Apr 2022 04:47:02 PM CST
################################################

#!/bin/sh 

set -eou pipefail

# script=/home/user/data2/lit/project/6mA/bin/binned.sh

# time ( bash $script \
# -s $(for i in `seq 1 2`;do printf "PC10_LY_MCC-M$i|";done) \
# -o /home/user/data/lit/project/6mA/feature/danpos/normed/binned/PC10_LY_MCC-M \
# -fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
# -ld log/danpos/normed/binned/PC10_LY_MCC-M ) 


# script=/home/user/data2/lit/project/6mA/bin/binned_level.sh

# time ( bash $script \
# -s $(for i in `seq 1 2`;do printf "PC10_LY_MCC-M$i|";done) \
# -o /home/user/data/lit/project/6mA/feature/danpos/normed/test \
# -mA "/home/user/data2/lit/project/6mA/mA_list/sig_gt_2_fc0_cov10_for_lt_wt1.sorted.bed|/home/user/data2/lit/project/6mA/mA_list/sig_gt_2_fc0_cov10_for_lt_wt2.sorted.bed|/home/user/data2/lit/project/6mA/mA_list/sig_gt_2_fc0_cov10_for_lt_wt1-wt2.sorted.bed" \
# -fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
# -fa /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa \
# -ld /home/user/data2/lit/project/6mA/log/danpos/normed/test )

script=/home/user/data2/lit/project/6mA/bin/binned_level.sh

time ( bash $script \
-s $(for i in `seq 1 2`;do printf "PC10_LY_MCC-M$i|";done) \
-o /home/user/data/lit/project/6mA/feature/danpos/normed/AT_rich_list \
-mA "/home/user/data2/lit/project/6mA/mA_list/sig_gt2_fc0_wt1-wt2_merge_correct-strand-1.sorted.bed|/home/user/data2/lit/project/6mA/mA_list/sig_gt2_fcnot0_wt1-wt2_merge_correct-strand-2.sorted.bed|/home/user/data2/lit/project/6mA/mA_list/merge_correct-strand-3.sorted.bed" \
-fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
-fa /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa \
-ld /home/user/data2/lit/project/6mA/log/danpos/normed/AT_rich_list \
-bd /home/user/data/lit/project/6mA/feature/danpos/normed/binned/PC10_LY_MCC-M )


# nohup bash /home/user/data2/lit/project/6mA/bin/run.binned_level.sh &