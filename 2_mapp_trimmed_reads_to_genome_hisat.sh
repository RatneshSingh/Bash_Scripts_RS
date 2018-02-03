function usage()
{
    echo "$0 Read_folder  Output_folder[2_trimmed_reads_mapped]"
    exit;
    
}


tread_fold=$1
map_fold=$2

[[ ! -z $tread_fold ]] || usage ;

: ${map_fold:=2_trimmed_reads_mapped_hisat}

wdir=$map_fold
data_dir="/home/ratnesh.singh/Turf/0_final_files"
genomeFasta="$data_dir/Diamond_bionanoscaffolds_LG_Allmaps.fasta"
read_dir=$(realpath $tread_fold)
#GFF_file=/mnt/nfs1/ratnesh.singh/Turf/0_final_files/Diamond_bionanoscaffolds_LG_Allmaps.fasta.evm.gffread.gtf
GFF_file=/home/ratnesh.singh/Turf/0_final_files/diamond_jamg.gene_structures_post_PASA_updates.63725.transcript.gtf 






mkdir $map_fold
mapfold=$(realpath $map_fold)
cd $map_fold


for files in $(ls ../$tread_fold/*R1*.P.qtrim.gz); do 
file_R1=$(realpath $files)
file_R2=$(realpath ${files/R1/R2})
readFasta=" $file_R1 $file_R2 "
reads=" $file_R1 $file_R2 "
out_prefix=$(basename $files)
out_prefix=${out_prefix%%R1*}
out_prefix=HISAT2_map_${out_prefix}_R1R2_P

reads_commnd="  --readFilesCommand zcat --readFilesIn $reads "




##############################################################################################


outTmpDir="/home/ratnesh.singh/local_storage/${out_prefix}_HISAT2tmp"
rm -r $outTmpDir

#genomeIndex=$wdir/$(basename $genomeFasta).star.index
genomeIndex=${genomeFasta}.1.ht2

### create genome index for STAR if does not exists.
if [[ -e "$genomeIndex" ]];then
    echo "Index file for $genomeFasta exists. Skipping index generation step"
else
    
    hisat2-build -p 40 --quiet --exon $GFF_file    $genomeFasta  $genomeFasta

fi
    

### Map transcripts to genome using hisat2
#### hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
mkdir $outTmpDir
nice hisat2   --max-intronlen 70000 --dta -p 46 -mm  -x $genomeFasta -m1 $file_R1 -m2 $file_R2 -S $outTmpDir/${out_prefix}.hisat.sam  && samtools sort -O bam  -@ 6 -o $outTmpDir/${out_prefix}.hisat.sorted.bam ${out_prefix}.hisat.sam && mv $outTmpDir/${out_prefix}.hisat.sorted.bam $map_fold/



done



