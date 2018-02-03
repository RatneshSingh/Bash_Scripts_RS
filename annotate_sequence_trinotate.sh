usage="$0 -t transcript 
Options...
-r Do not RNAMer stage. It giver error when no tRNA is found.
-n Do not treat transcript as Trinity output.

"


while getopts ":t:rn" opt; do 
case $opt in
 t)
   transcript=$OPTARG
   ;;
 r)
   rnammer=1
  ;;
 n)
  no_trinity=1
 ;;
 \?)
   echo "Invalid options: -$OPTARG" >&2
   exit 1
   ;;
 
  :)
   echo "Option -$OPTARG requires an argument." >&2
   exit 1

esac
done


if [ $OPTIND -eq 1 ]; then echo "Script needs options";echo; echo "$usage"; exit 1; fi
shift $((OPTIND-1))
[[ $# -eq 0 ]] || echo "$# non-option arguments were used."



if [[ ! -s $transcript ]]; then echo "$transcript not found" ; exit 1; fi 
#exit;
###### HOME folders for used programs #############
TRINITY_HOME="/usr/local/trinity"
TRINOTATE_HOME="/usr/local/trinotate"
DATABASE_HOME="/mnt/nfs1/Databases/trinotate_databases/TRINOTATE_V3"
temp_folder="/home/${USER}/local_storage/tmp_folders"
cpu=54
cross_sample_norm="TMM"   ### cross sample normalization methods e.g.  TMM, UpperQuartile, none

transcript=$(realpath $transcript)
anno_fold=${transcript%/*}
transcript_name=$(basename $transcript)
transcript_folder=$(dirname $transcript)
transcript=$transcript_folder/$transcript_name
anno_dir=$transcript_folder/${transcript_name}_Trinotate_annotation


###################### choose conf file for Trinotate
## choose this one no rnammer run.
trino_conf=$TRINOTATE_HOME/auto/conf_diamond_no_rnammer.txt
## choose this one if you want with rnammer. It causes problem when there is no rna found.
#trino_conf=$TRINOTATE_HOME/auto/conf_diamond.txt


### 
curdir=$PWD
mkdir -p $anno_dir
cd $anno_dir
trnsname=$transcript_name
[[ -e $transcript_name ]] || ( ln -r -s ${transcript} $trnsname && ln -s -f -r  ${transcript} Trinity.fasta )



if [[ ! -f $trnsname.gene_trans_map || ! -s $trnsname.gene_trans_map ]];then
 #$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl $trnsname >  $trnsname.gene_trans_map
 echo "Not Running as Trinity. If it is Trinity generated transcript. uncomment the trns generationstep using trinity"
fi

 if [[ ! -s $trnsname.gene_trans_map ]]; then
	grep ">" $trnsname|sed 's|>||g'|awk '{printf ("%s\t%s\n",$1,$1)}' > $trnsname.gene_trans_map
 fi


### create a symlink to gene length file.
if [[ ! -f $trnsname.gene_iso.seq_lengths ]];then
  [[ -e ${transcript}.gene_iso.seq_lengths ]] && ln -r -s ${transcript}.gene_iso.seq_lengths  ${trnsname}.gene_iso.seq_lengths
  [[ ! -e ${trnsname}.gene_iso.seq_lengths ]] && perl  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence;} foreach(keys %sequence){print "$_\t".length($sequence{$_})."\n"} ' < ${transcript} > ${trnsname}.gene_iso.seq_lengths

fi
############################
### to speed up run following in parallel
# progs
TRANSDECODER_DIR=/usr/local/TransDecoder
BLASTX_PROG=blastx
BLASTP_PROG=blastp
DIAMOND_PROG=diamond
SIGNALP_PROG=/usr/local/signalp/signalp
TMHMM_PROG=/usr/local/tmhmm-2.0c/bin/tmhmm
RNAMMER_TRANS_PROG=$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl
RNAMMER=/usr/local/rnammer-1.2/rnammer
HMMSCAN_PROG=hmmscan

# dbs
SWISSPROT_PEP=/mnt/nfs1/Databases/trinotate_databases/TRINOTATE_V3/uniprot_sprot.pep
PFAM_DB=/mnt/nfs1/Databases/trinotate_databases/TRINOTATE_V3/Pfam-A.hmm


####################
#  BioIfx computes  ** no need to edit the commands below **
####################
#[[ -e trinotate_primary_steps_transdecoder_inparallel.txt ]] && rm trinotate_primary_steps_transdecoder_inparallel.txt

section_name=TRANSDECODER_LONGORF
[[ -e chkpt.$section_name.ok ]] || $TRANSDECODER_DIR/TransDecoder.LongOrfs -t $trnsname && $TRANSDECODER_DIR/TransDecoder.Predict -t $trnsname --cpu $cpu && touch "chkpt.$section_name.ok"

section_name=TRANSDECODER_PREDICT
[[ -e chkpt.$section_name.ok ]] || $TRANSDECODER_DIR/TransDecoder.Predict -t $trnsname --cpu $cpu && touch "chkpt.$section_name.ok"


### will need predcited ORF and peptides for further runs. so run it first.
[[ -e chkpt.TRANSDECODER_PREDICT.ok && -e chkpt.TRANSDECODER_LONGORF.ok ]] || (echo "Need TRANSDECODER run before other programs could be run" &&  exit)


[[ -e trinotate_primary_steps_inparallel.tx ]] && rm trinotate_primary_steps_inparallel.tx

section_name=BLASTX_SPROT_TRANS
[[ -e chkpt.$section_name.ok ]] || echo "$DIAMOND_PROG  blastx -d $SWISSPROT_PEP -q $trnsname -p $cpu -k 1  -e 0.00001 -a swissprot.blastx --sensitive && diamond view -a swissprot.blastx.daa -o swissprot.blastx.outfmt6 && touch "chkpt.$section_name.ok""|tee  trinotate_primary_steps_inparallel.txt

section_name=BLASTX_SPROT_PEP
[[ -e chkpt.$section_name.ok ]] || echo "$DIAMOND_PROG  blastp -d $SWISSPROT_PEP -q $trnsname.transdecoder.pep  -p $cpu -k 1  -e 0.00001 -a swissprot.blastp --sensitive ; diamond view -a swissprot.blastp.daa -o swissprot.blastp.outfmt6 && touch "chkpt.$section_name.ok""|tee -a trinotate_primary_steps_inparallel.txt

section_name=PFAM
[[ -e chkpt.$section_name.ok ]] || echo "$HMMSCAN_PROG --cpu $cpu --domtblout TrinotatePFAM.out $PFAM_DB $trnsname.transdecoder.pep  > pfam.log && touch "chkpt.$section_name.ok""|tee -a trinotate_primary_steps_inparallel.txt

section_name=SIGNALP
[[ -e chkpt.$section_name.ok ]] || echo "$SIGNALP_PROG -f short -n signalp.out $trnsname.transdecoder.pep > sigP.log && touch "chkpt.$section_name.ok""|tee -a trinotate_primary_steps_inparallel.txt

section_name=TMHMM
[[ -e chkpt.$section_name.ok ]] || echo "$TMHMM_PROG --short < $trnsname.transdecoder.pep > tmhmm.out && touch "chkpt.$section_name.ok" "|tee -a trinotate_primary_steps_inparallel.txt

section_name=RNAMMER
[[ -e chkpt.$section_name.ok ]] || echo "$RNAMMER_TRANS_PROG --transcriptome $trnsname --path_to_rnammer $RNAMMER  && touch "chkpt.$section_name.ok" "|tee -a trinotate_primary_steps_inparallel.txt
# generates file: $trnsname.rnammer.gff



parallel -j 8 < trinotate_primary_steps_inparallel.txt

###### Now run the pipeline with checkpoints and results in place to load.
if [[ ! -s $trnsname.rnammer.gff ]];then
 trino_conf=$TRINOTATE_HOME/auto/conf_diamond_no_rnammer.txt
else
   ## choose this one if you want with rnammer

  trino_conf=$TRINOTATE_HOME/auto/conf_diamond.txt
fi



$TRINOTATE_HOME/auto/autoTrinotate.pl --Trinotate_sqlite $DATABASE_HOME/Trinotate.sqlite  --transcripts $trnsname --gene_to_trans_map $trnsname.gene_trans_map  --conf $trino_conf --CPU $cpu

$TRINOTATE_HOME/util/Trinotate_get_feature_name_encoding_attributes.pl Trinotate.xls > annot_feature_map.txt

$TRINOTATE_HOME/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -T  --include_ancestral_terms > go_annotations_isoforms.txt
$TRINOTATE_HOME/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -G  --include_ancestral_terms > go_annotations_genes.txt

$TRINITY_HOME/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl go_annotations_genes.txt annot_feature_map.txt > wAnnot.go_annotations_genes.txt
$TRINITY_HOME/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl go_annotations_isoforms.txt annot_feature_map.txt > wAnnot.go_annotations_isoforms.txt

cat go_annotations_genes.txt go_annotations_isoforms.txt > go_annotations.txt
cat wAnnot.go_annotations_genes.txt wAnnot.go_annotations_isoforms.txt > wAnnot.go_annotations.txt


[[ -e wAnnot.${trnsname}.gene_iso.seq_lengths ]] || $TRINITY_HOME/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl ${trnsname}.gene_iso.seq_lengths annot_feature_map.txt > wAnnot.${trnsname}.gene_iso.seq_lengths




cd $curdir
