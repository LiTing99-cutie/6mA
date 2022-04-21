#!/bin/sh

set -eou pipefail

<<'!'
    script=/home/user/data2/lit/project/6mA/reproduce/bin/6mAscope.sh

    bash ${script} \
    /bam/Algae.subreads.bam \
    /home/user/data2/lit/project/6mA/reproduce/output/Pacbio \
    /home/user/data2/lit/project/6mA/reproduce/log &

    bash ${script} \
    /bam/Tetrahymena.subreads.bam \
    /home/user/data2/lit/project/6mA/reproduce/output/Pacbio \
    /home/user/data2/lit/project/6mA/reproduce/log &

    # not work
    # wait
!

<<'!'
script=/home/user/data2/lit/project/6mA/reproduce/bin/6mAscope.ipd.sh

bash ${script} \
/bam/Algae.subreads.bam \
/home/user/data2/lit/project/6mA/reproduce/output/Pacbio \
/home/user/data2/lit/project/6mA/reproduce/log &

# bash ${script} \
# /bam/Tetrahymena.subreads.bam \
# /home/user/data2/lit/project/6mA/reproduce/output/Pacbio \
# /home/user/data2/lit/project/6mA/reproduce/log &
!


script=bin/6mAscope.sh

time ( bash ${script} bam/P05DY20245308-1_r64118_20200818_093653_3_G01.subreads.bam ) >SMRT-1.6mASCOPE.log 2>&1

# nohup bash /home/user/data2/lit/project/6mA/reproduce/bin/run.6mAscope.sh &