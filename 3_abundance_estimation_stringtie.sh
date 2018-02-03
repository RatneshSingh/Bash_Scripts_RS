function exit_on_error(){

echo "$@"

echo "Usage: sh $0         bam_folder       out_folder    ref_gtf  trnscript_label[STRG]   pattern_to_grep_for_file[sorted.bam]"

exit;

}


bam_fold=$1
out_fold=$2
ref_gtf=$3
label=$4
patt=$5

: ${ref_gtf:=~/Turf/0_final_files/diamond_jamg.gene_structures_post_PASA_updates.63725.transcript.gtf}
: ${label:="STRG"}
wdir=$PWD
curdir=$(dirname $wdir)

[[ ! -z "$bam_fold" ]] ||  exit_on_error "Provide path for bam folder or file"

: ${out_fold:=3_abundance_estimation_stringtie}
: ${patt:=sorted.bam}

[[ -d ${out_fold} ]] || mkdir -p $out_fold

declare -a file_list

if [ -f "${bam_fold}" ]; then
  echo "$bam_fold is a file."
  file_list=($(realpath $file))
elif [ -d "${bam_fold}" ]; then
  echo "$bam_fold is a folder. looking for $patt containing bam files"
  file_list=($(ls -d1  $(realpath $bam_fold)/*.* | grep $patt | grep "bam$"))
fi


echo "Working on files ${file_list[@]}"


for file in "${file_list[@]}"; do  
        filename=$(basename $file)	
	#label=${filename:0:20}
	[[ -s $wdir/$out_fold/${filename}.stringtie.abund.out ]] && { >&2 echo "outfile $wdir/$out_fold/${filename}.stringtie.abund.out exists. Delete before running this script"; continue; }
        echo "Running stringtie on $file using parallel"
        #echo "stringtie $file -G ${ref_gtf} -l $label -p 40 -A $wdir/$out_fold/${filename}.stringtie.abund.out -x  -e -b $wdir/$out_fold/${filename}.forBallgown_outdir -o $wdir/$out_fold/${filename}.stringtie.assembly.gtf" 
        echo "stringtie $file -G ${ref_gtf} -p 40 -A $wdir/$out_fold/${filename}.stringtie.abund.out -x  -e -b $wdir/$out_fold/${filename}.forBallgown_outdir -o $wdir/$out_fold/${filename}.stringtie.assembly.gtf" 

done | parallel -j 20



