#!/bin/sh 

################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/run/run.binned.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 18 Apr 2022 04:50:28 PM CST
################################################


set -eou pipefail

script=/home/user/data2/lit/project/6mA/bin/binned.update.sh

logDir=/home/user/data2/lit/project/6mA/reproduce/log

<<'!'
time (
bash ${script} \
C.reinhardtii \
/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/C.reinhardtii \
/home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
/home/user/data/lit/project/6mA/reproduce/iNPS_output/C.reinhardtii_Gathering.like_bed \
TRUE \
10 
) > $logDir/run.binned.C.reinhardtii.log 2>&1 &

time (
bash ${script} \
T.thermophila \
/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/T.thermophila \
/home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
/home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed \
TRUE \
10
) > $logDir/run.binned.T.thermophila.log 2>&1 &

for sample in PC10_LY_MCC-M1 PC10_LY_MCC-M2;do
time (
bash ${script} \
$sample \
/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/$sample \
/home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
/home/user/data2/lit/project/6mA/danpos/normed/$sample/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls \
TRUE \
10
) > $logDir/run.binned.$sample.log 2>&1 &
done
!

bin_number=10
# script
    time (
    bash ${script} \
    C.reinhardtii \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/C.reinhardtii \
    /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/C.reinhardtii_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.C.reinhardtii.${bin_number}bin.log 2>&1 &

    time (
    bash ${script} \
    T.thermophila \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/T.thermophila \
    /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.T.thermophila.${bin_number}bin.log 2>&1 &

    for sample in PC10_LY_MCC-M1 PC10_LY_MCC-M2;do
    time (
    bash ${script} \
    $sample \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/$sample \
    /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
    /home/user/data2/lit/project/6mA/danpos/normed/$sample/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.$sample.${bin_number}bin.log 2>&1 &
    done

# ---

<<'!'
bin_number=8
# script
    time (
    bash ${script} \
    C.reinhardtii \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/C.reinhardtii \
    /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/C.reinhardtii_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.C.reinhardtii.${bin_number}bin.log 2>&1 &

    time (
    bash ${script} \
    T.thermophila \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/T.thermophila \
    /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.T.thermophila.${bin_number}bin.log 2>&1 &

    for sample in PC10_LY_MCC-M1 PC10_LY_MCC-M2;do
    time (
    bash ${script} \
    $sample \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/$sample \
    /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
    /home/user/data2/lit/project/6mA/danpos/normed/$sample/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.$sample.${bin_number}bin.log 2>&1 &
    done

# ---

bin_number=6
# script
    time (
    bash ${script} \
    C.reinhardtii \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/C.reinhardtii \
    /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/C.reinhardtii_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.C.reinhardtii.${bin_number}bin.log 2>&1 &

    time (
    bash ${script} \
    T.thermophila \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/T.thermophila \
    /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.T.thermophila.${bin_number}bin.log 2>&1 &

    for sample in PC10_LY_MCC-M1 PC10_LY_MCC-M2;do
    time (
    bash ${script} \
    $sample \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/$sample \
    /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
    /home/user/data2/lit/project/6mA/danpos/normed/$sample/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.$sample.${bin_number}bin.log 2>&1 &
    done

# ---

bin_number=4
# script
    time (
    bash ${script} \
    C.reinhardtii \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/C.reinhardtii \
    /home/user/data/lit/database/public/genome/C.reinhardtii/Creinhardtii_281_v5.0.fa.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/C.reinhardtii_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.C.reinhardtii.${bin_number}bin.log 2>&1 &

    time (
    bash ${script} \
    T.thermophila \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/T.thermophila \
    /home/user/data/lit/database/public/genome/T.thermophila/T_thermophila_June2014_assembly.fasta.sorted.size \
    /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.T.thermophila.${bin_number}bin.log 2>&1 &

    for sample in PC10_LY_MCC-M1 PC10_LY_MCC-M2;do
    time (
    bash ${script} \
    $sample \
    /home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/${bin_number}bin/$sample \
    /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size \
    /home/user/data2/lit/project/6mA/danpos/normed/$sample/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls \
    TRUE \
    ${bin_number}
    ) > $logDir/run.binned.$sample.${bin_number}bin.log 2>&1 &
    done

# ---
!
wait