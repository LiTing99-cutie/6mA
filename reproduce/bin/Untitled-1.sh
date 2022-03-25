
pushd FileFromWindows

unzip download.20220307.210139.zip

cp Phytozome/PhytozomeV10/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene.gff3.gz $gtf_dir
cp Phytozome/PhytozomeV10/Creinhardtii/assembly/Creinhardtii_281_v5.0.fa.gz $genome_data_dir_c

rm -rf download.20220307.210139.zip Phytozome/

popd

# FileFromWindows/Download_78086_File_Manifest.csv are metadata of Creinhardtii data