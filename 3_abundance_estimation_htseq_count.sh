function exit_on_error(){

echo "$@"

echo "Usage: $0 bam_folder out_folder pattern_to_grep_for_file[sorted.bam]"

exit;

}


bam_fold=$1
out_fold=$2
patt=$3

#ref_gtf=~/Turf/0_final_files/diamond_jamg.gene_structures_post_PASA_updates.63725.transcript.gtf
ref_gtf=~/Turf/0_final_files/diamond_AllRNASeq_stringtie_assembled_merged.gtf
wdir=$PWD
curdir=$(dirname $wdir)

[[ ! -z "$bam_fold" ]] ||  exit_on_error "Provide path for bam folder or file"

: ${out_fold:=3_abundance_estimation_htseq_count}
: ${patt:=sorted.bam}

mkdir -p $wdir/$out_fold

declare -a file_list

if [ -f "${bam_fold}" ]; then
  echo "$bam_fold is a file."
  file_list=($(realpath $file))
elif [ -d "${bam_fold}" ]; then
  echo "$bam_fold is a folder. looking for bam files matching $patt in names"
  file_list=($(ls -d1  $(realpath $bam_fold)/*.* | grep $patt | grep "bam$"))
fi


#echo "Working on files $file_list"


for file in "${file_list[@]}"; do  
        filename=$(basename $file)
        [[ -s $wdir/$out_fold/${filename}.htseqcount_nonstranded.abund.out ]] && { >&2 echo "Count file exists. Delete it to run htseq on it again.  Skipping $filename"; continue; }	
        echo "Running htseq for $file in detached screen in background."
	screen -m -d -S run_htseq_count_detached bash -c "htseq-count -r pos -s no -f bam  $file ${ref_gtf} >  $wdir/$out_fold/${filename}.htseqcount_nonstranded.abund.out"  ; 

done



