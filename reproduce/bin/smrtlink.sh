#!/bin/sh

set -eou pipefail
# conda create -n pbalign -c bioconda pbalign

. "/home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh"
conda activate pbalign

subreadBam=$1
basename=`basename $subreadBam`
subreadBamName=${basename%%.*}
ref=$2
ccs_fasta=$3
output_dir=$4/${subreadBamName}
log_dir=$5

[ -d ${output_dir} ] || mkdir ${output_dir}
[ -d ${log_dir} ] || mkdir ${log_dir}

# minimap2 2.15-r905
# minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]

# pbalign 0.4.1
# positional arguments:
#   inputFileName         SubreadSet or unaligned .bam
#   referencePath         Reference DataSet or FASTA file
#   outputFileName        Alignment results dataset
# Input can be a
# fasta, pls.h5, bas.h5 or ccs.h5 file or a fofn (file of file names). Output
# can be in SAM or BAM format. If output is BAM format, aligner can only be
# blasr and QVs will be loaded automatically. NOTE that pbalign no longer
# supports CMP.H5 in 3.0.
# The input file format can only be FASTA,PLS_H5,PLX_H5,BAS_H5,BAX_H5,FOFN,CCS_H5,BAM,XML. (from log)

# ipdSummary 3.0
# ipdSummary aligned.bam --reference ref.fasta m6A,m4C --gff basemods.gff \
# --csv_h5 kinetics.h5

time ( pbalign --nproc 50 \
--tmpDir /home/user/data2/lit/project/6mA/reproduce/tmp \
$subreadBam \
${ccs_fasta} \
${output_dir}/${subreadBamName}.pbalign.ccs.bam ) \
> ${log_dir}/pbalign.${subreadBamName}.ccs.log 2>&1 

time ( ipdSummary --numWorkers 50 \
${output_dir}/${subreadBamName}.pbalign.ccs.bam \
--reference ${ccs_fasta} \
--identify m6A,m4C \
--methylFraction \
--gff ${output_dir}/ccs/${subreadBamName}.gff \
--csv ${output_dir}/ccs/${subreadBamName}.csv \
--bigwig ${output_dir}/ccs/${subreadBamName}.bigwig ) \
> ${log_dir}/ipdSummary.${subreadBamName}.ccs.log 2>&1 

# time ( pbalign --nproc 50 \
# --tmpDir /home/user/data2/lit/project/6mA/reproduce/tmp \
# $subreadBam \
# $ref \
# ${output_dir}/${subreadBamName}.pbalign.ref.bam ) \
# > ${log_dir}/pbalign.${subreadBamName}.ref.log 2>&1 

# time ( ipdSummary --numWorkers 50 \
# ${output_dir}/${subreadBamName}.pbalign.ref.bam \
# --reference ${ref} \
# --identify m6A,m4C \
# --methylFraction \
# --gff ${output_dir}/ref/${subreadBamName}.gff \
# --csv ${output_dir}/ref/${subreadBamName}.csv \
# --bigwig ${output_dir}/ref/${subreadBamName}.bigwig ) \
# > ${log_dir}/ipdSummary.${subreadBamName}.ref.log 2>&1 

