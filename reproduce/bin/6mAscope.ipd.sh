# #!/bin/sh

subreadBam=$1
basename=`basename $subreadBam`
subreadBamName=${basename%%.*}
output_dir=$2/${subreadBamName}
log_dir=$3

[ -d ${output_dir} ] || mkdir -p ${output_dir}
[ -d ${log_dir} ] || mkdir -p ${log_dir}

pushd $output_dir

time ( 6mASCOPE ipd -r $subreadBam -f $subreadBamName.ccs.fasta -o $subreadBamName.Ipd.out -p 40 ) \
>${log_dir}/ipd.${subreadBamName}.log 2>&1

popd
