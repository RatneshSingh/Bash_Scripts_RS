repfile=$1
reffile=$2
len=$3
span=$4
pad=$5
mfe=$6
maln=$7


bamfold=$8


#len=25  ### length of unmasked region in the sequence.
#span=10 ### Read should span this many bases accross the junction on both sides to be counted as present in accession.
#pad=5   ### pad the ends to merge hits withinh a window of this size.
###################################################################################
if [  -z $repfile ];then
    echo;echo;echo "Usage: script repeat_file reference_file unmasked_length span_junc pad_locus maxFromEnd MinAln bamfolder"
    echo;echo;exit
fi

### mask internal region of each mite leaving 25bp at the ends.

perl -w -e '
$/="\n>"; 
my$file=$ARGV[0];
my$mask_len=$ARGV[1];
open REP,"$ARGV[0]";
while(<REP>){
my($seq,$nseq,$header,@sequence);
#print "** original: $_\n";
($header,@sequence)=split /\n/,$_;
$seq=join("",@sequence);
$seq=~s/\s*|>|[^ATGCNatgcn]//g;
$header=~s/>//g;
my$num_Ns= (length($seq)- 2*$mask_len)  > 99 ? length($seq) -  2 * $mask_len : 100;

substr($seq,$mask_len -1 ,-($mask_len), "N" x $num_Ns); 
print ">$header\n$seq\n";   
}' $repfile $len > ${repfile%.*}.$len.fa



#### blast masked repeat sequences to reference genome
#if [ -e "$reffile" ]; then
 
#### make blastdb if it does not exists
# if [ ! -e "$reffile.nhr" ]; then 
#    makeblastdb -in $reffile -dbtype nucl
#fi
 ### run blastn if file does not exists
 #if [ ! -e ${repfile%.*}.$len.vs.${reffile%.*}.blastn.table ]; then
 #  blastn -task blastn  -query ${repfile%.*}.$len.fa -db $reffile -out ${repfile%.*}.$len.vs.${reffile%.*}.blastn.table -outfmt 6  -num_threads 20
 #fi


 if [ ! -e "${repfile%.*}.$len.fa.nhr" ]; then
     echo "Creating Blast database for ${repfile%.*}.$len.fa"
     makeblastdb -in ${repfile%.*}.$len.fa -dbtype nucl
  else
     echo "Formatted database for $repfile exists. Skipping makeblastdb run"
 fi
 
 if [ ! -e ${reffile%.*}.vs.${repfile%.*}.$len.blastn.table ]; then
     echo "Running blast for"
     echo "blastn -task blastn  -query $reffile -db ${repfile%.*}.$len.fa -out ${reffile%.*}.vs.${repfile%.*}.$len.blastn.table -outfmt '6 std qlen slen' -culling_limit 1 -xdrop_gap 50 -max_target_seqs 10000000"
     blastn -task blastn  -query $reffile -db ${repfile%.*}.$len.fa -out ${reffile%.*}.vs.${repfile%.*}.$len.blastn.table -outfmt '6 std qlen slen' -culling_limit 1 -xdrop_gap $len -max_target_seqs 10000000 -num_threads 70
     blastn -task blastn  -query $reffile -db ${repfile%.*}.$len.fa -out ${reffile%.*}.vs.${repfile%.*}.$len.blastn -culling_limit 1 -xdrop_gap $len -max_target_seqs 10000000 -num_threads 70 &
 else
     echo "Blast result already exists. Skipping blastn"
 fi

### assign blast result for later use
blast_table=${reffile%.*}.vs.${repfile%.*}.$len.blastn.table


#### create a non-redundant table for all the hits
#if [ -e "$blast_table" ]; then

#echo "Table redundant" 
#perl -e 'while(<>){}' 
#  perl ~/Scripts/Perl/Blast_Table_Sumup_Redundant_hits_V3.0.pl -b ${repfile%.*}.$len.vs.${reffile%.*}.blastn.table  -o 

#fi 

#### create n=bam file list
#printf "### file names:" > ${repfile%.*}.$len.vs.${reffile%.*}.blastn.table.subject.summary.table.readcount.table
#for file in $bamfold/*.bam;do
#echo $file;
#printf "\t$file" >> ${repfile%.*}.$len.vs.${reffile%.*}.blastn.table.subject.summary.table.readcount.table
#done > bam_files.list



#### convert to bed format and get the number of reads passingthrough the junction (+- 10 bp)
## create the commandlist to run using parallel
for file in $bamfold/*.bam; do
    bamname=${file##*/};
    if [ ! -e ${file/bam/bed} ]; then
       echo bedtools bamtobed -i $file \> ${file/bam/bed};
    fi
       
done > bamtobed_commands.list 

#parallel -j 10 -a bamtobed_commands.list 


perl 2_extract_repeat_junctions.pl  -bls $blast_table -sae $len -spj $span -pad $pad -mfe $mfe -maln $maln -bff $bamfold/*.bed






