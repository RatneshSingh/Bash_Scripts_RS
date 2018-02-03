function usage()
{
    echo "$0 Read_folder  Output_folder[2_trimmed_reads_mapped]"
    exit;
    
}


tread_fold=$1
map_fold=$2

[[ ! -z $tread_fold ]] || usage ;

: ${map_fold:=2_trimmed_reads_mapped}

wdir=$map_fold
data_dir="/home/ratnesh.singh/Turf/0_final_files"
genomeFasta="$data_dir/Diamond_bionanoscaffolds_LG_Allmaps.fasta"
read_dir=$(realpath $tread_fold)
#GFF_file=/mnt/nfs1/ratnesh.singh/Turf/0_final_files/Diamond_bionanoscaffolds_LG_Allmaps.fasta.evm.gffread.gtf
GFF_file=/home/ratnesh.singh/Turf/0_final_files/diamond_jamg.gene_structures_post_PASA_updates.63725.transcript.gtf 






mkdir $map_fold
cd $map_fold


for files in $(ls ../$tread_fold/*R1*.P.qtrim.gz); do 
file_R1=$(realpath $files)
file_R2=$(realpath ${files/R1/R2})
readFasta=" $file_R1 $file_R2 "
reads=" $file_R1 $file_R2 "
out_prefix=$(basename $files)
out_prefix=${out_prefix%%R1*}
out_prefix=STAR_map_${out_prefix}_R1R2_P

reads_commnd="  --readFilesCommand zcat --readFilesIn $reads "




##############################################################################################


outTmpDir="/home/ratnesh.singh/local_storage/${out_prefix}_STARtmp"
outTmpDir="/dev/shm/${out_prefix}_STARtmp"
outfile=${out_prefix}Aligned.out.bam
[[ -s $outfile ]] && echo "$map_fold/$(basename $outfile) exists. Skipping..."
[[ -s $outfile ]] && continue
[[ -d $outTmpDir ]] && echo "$outTmpDir exists. Skipping..."
[[ -d $outTmpDir ]] && continue



#rm -r $outTmpDir


#continue


#genomeIndex=$wdir/$(basename $genomeFasta).star.index
genomeIndex=${genomeFasta}.star.index

### create genome index for STAR if does not exists.
if [[ -e "$genomeIndex/SAindex" ]];then
    echo "Index file for $genomeFasta exists. Skipping index generation step"
else
    mkdir -p $genomeIndex    
    STAR \
    --runThreadN 40 \
    --runMode genomeGenerate \
    --genomeDir $genomeIndex \
    --sjdbGTFfile $GFF_file \
    --genomeFastaFiles $genomeFasta 

fi
    

### Map transcripts to genome using STAR.
#nice STARlong \
#--runThreadN 40 \
#--genomeDir $genomeIndex \
#--genomeLoad Remove

echo "will run STAR on $reads_commnd"

#continue
#exit

nice STARlong \
--runThreadN 40 \
--genomeDir $genomeIndex \
--runMode alignReads \
--outFilterType BySJout \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--outFilterMultimapNmax 20 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1\
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--quantMode TranscriptomeSAM GeneCounts \
--sjdbGTFfile $GFF_file \
--outFileNamePrefix $out_prefix \
--outTmpDir $outTmpDir \
--twopassMode Basic \
$reads_commnd

#--outSAMtype BAM SortedByCoordinate \
#--twopassMode Basic \
#--sjdbFileChrStartEnd ${out_prefix}_STARpass1/SJ.out.tab \
#--sjdbGTFtagExonParentTranscript Parent \


[[ -s $outfile ]] || exit
#cp -r  $outTmpDir/*   $map_fold/ && rm -r $outTmpDir

done



