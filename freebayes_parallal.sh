curdir=$PWD
wdir="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations"
restitution_folder=/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/11_genome_restitution_on_LAP9202Ref

refname="LAP.US56.combined"  ### will be used to name folder and files.
transcript="$wdir/2_trimmed_reads_merged/${refname}/LAP.US56.9202.combined.fasta"
estmeth="RSEM"
TRINITY_HOME="/usr/local/trinity"

read_folder="2_trimmed_reads"
#read_folder="1_raw_reads"

trinity="${TRINITY_HOME}/Trinity"
temp_folder="/home/ratnesh.singh/local_storage/vcf_folders"
local_wdir="/home/ratnesh.singh/local_storage"
matrix_fold="$wdir/7.1_Expression_Matrix"

ploidy=12
grp="All"
min_coverage=2
ref=$transcript
bam=$bam_list
vcfout="$restitution_folder/${grp}.All_sorted.rg.freebayes_2ndRun.vcf"
use_best_n_alleles=8
min_coverage=3
F=0.1
num_jobs=64

[[ -e $ref.fai ]] || /usr/local/bin/samtools  faidx  $ref

if [[ -e $ref.fai ]]; then
[[ -e ${ref}.regions ]] || fasta_generate_regions.py $ref.fai 100000 > ${ref}.regions
fi


num_lines=$(wc -l ${ref}.regions)
per_job=$(( $num_lines/$num_jobs + 1 ))

for job in {0..${num_jobs}; do
    start=$(($num_jobs * $per_job + 1))
    end=$(($start+$perjob))
    region_out=${ref}.regions.${job}_${per_job}
    
    [[ -e $region_out ]] || tail -n +$start ${ref}.regions|head -n $per_job > $region_out
    [[ -e $region_out.vcfout.working || -e $region_out.vcfout.ok ]] || (touch ${ref}.regions.${job}_${per_job}.vcfout.working && freebayes --targets $region_out --fasta-reference $ref --ploidy $ploidy --use-best-n-alleles ${use_best_n_alleles} --pooled-continuous --min-coverage $min_coverage -F $F --vcf ${ref}.regions.${job}_${per_job}.vcfout $bam && touch ${ref}.regions.${job}_${per_job}.vcfout.ok)
    
    rm touch ${ref}.regions.${job}_${per_job}.vcfout.working

done
    
    



freebayes --targets --fasta-reference $ref --ploidy $ploidy --use-best-n-alleles $use-best-n-alleles --pooled-continuous --min-coverage $min_coverage -F $F --vcf $vcfout $bam

&& echo "Done " && touch $restitution_folder/${grp}.All_sorted.rg.freebayes.vcf.ok