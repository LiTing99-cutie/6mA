################################################
#File Name: bin/run.filter_fuzziness_score_bin_level.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Mar 2022 04:10:15 PM CST
################################################

#!/bin/sh 
set -e

### wig2bed ###
# printf no line break
for sample in $(for i in `seq 1 2`;do echo "PC10_LY_MCC-M$i";done);do
wig2bed <  danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.wig > \
danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.wig.bed
done
### EO wig2bed ###

# annotation (single line) is not permitted in bash script
script=/home/user/data2/lit/project/6mA/bin/filter_fuzziness_score_bin_level.sh
time ( bash $script \
-s $(for i in `seq 1 2`;do printf "PC10_LY_MCC-M$i|";done) \
-o /home/user/data/lit/project/6mA/feature/danpos/normed/method_3 \
-mA "sig_gt_2_fc0_cov10_for_lt_wt1.bed|sig_gt_2_fc0_cov10_for_lt_wt2.bed|sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed" \
-c 100 \
-fs /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
-fa /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa \
-ld log/danpos/normed/method_3 ) 

# nohup bash bin/run.filter_fuzziness_score_bin_level.sh &


