################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/bin/run.smrtlink.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 24 Mar 2022 11:08:15 AM CST
################################################

#!/bin/sh 

ccs \
--min-rq 0.99 \
--min-passes 20 \

# For each DNA molecule, Circular Consensus reads (CCS reads) were generated using CCS software (v4.0.0) with 
# --min-rq 0.99 (i.e. –minPredictedAccuracy 0.99 in older CCS version).
# Only CCS reads with passes ≥ 20 were kept for the following analysis. 

# Subreads were hereafter mapped to the CCS read with pbalign (v0.4.1) with default parameters.

# IPD ratio and modification QV were calculated with ipdSummary (v2.4.1) from the mapping results.

# CCS reads were mapped to their reference genomes using minimap2 (v2.17-r941) (48) with –x sr, separately.