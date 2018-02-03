
function Cal_contamination(){
 local  sequencefile=$1
 local  blastdatabase=$2
 local  Resultfile=$3
 local  workdir=$4
 local extra_blast=$5

 ## create names for seq, db and outfiles
  local dbname=$(basename $blastdatabase)
  local seqname=$(basename $sequencefile)
  local outname=${seqname%.*}.vs.${dbname%.*}

  ## check if result file exists. If not, create new and add header.
  if [ ! -f $Resultfile ];then
    printf "Sequence_name\tTotal_seq_num\tTotal_seq_length\tClean_seq_num\tClean_seq_length\tContamin_seq_num\tContamin_seq_length\n" > $Resultfile
  fi

### if the database file id not formatted for blast, create blast db.
 if [ ! -f ${blastdatabase}.nhr ] && [ -z $no_format ] ; then
   makeblastdb -in $blastdatabase -dbtype nucl
 fi

  cd $workdir

  if [ -z $noblast ]; then
   echo "Running blast for $sequencefile against $blastdatabase"
   blastn -query $sequencefile -db $blastdatabase -out $workdir/$outname.blastn -num_threads 30 -evalue 1e-10 -xdrop_gap 3000 -xdrop_gap_final 5000 -culling_limit 2 -num_alignments 30000  #-outfmt '6 std qcovs stitle ' $extra_blast

  fi

  echo "Running blastparser for $workdir/$outname.blastn"
  perl $scriptfolder/Blast_parser_BIOPERL_V1.0.pl -i $workdir/$outname.blastn -o $workdir/$outname.blastn.p$perid.c$qcov.table -p $perid -c $qcov
  #perl $scriptfolder/

  echo "Getting list of sequence names showing hits to $dbname"
  cut -f1 $workdir/$outname.blastn.p$perid.c$qcov.table|sort|uniq > $workdir/$outname.blastn.p$perid.c$qcov.list
  echo "Removing contamination sequences from $seqname"
  perl $scriptfolder/Remove_selectedSequences-V1.3.pl -s $sequencefile -l $workdir/$outname.blastn.p$perid.c$qcov.list -m exact -d space -c 1

    ### count total bases
  perl $scriptfolder/Clean_FastaSeq_of_Numbers-V1.0.pl $sequencefile ${sequencefile}.clean
  totalnumseq=""
  totalnumseq=`grep -v ">" ${sequencefile}.clean|wc -lc|sed -r 's/\s+/\t/g'|sed -r 's/^\s+//g'`
  rm ${sequencefile}.clean

    ### count the clean bases
  perl $scriptfolder/Clean_FastaSeq_of_Numbers-V1.0.pl $sequencefile.include.fasta ${sequencefile}.include.fasta.clean
  mv ${sequencefile}.include.fasta.clean $workdir/$seqname.cleanedOf.$dbname.p$perid.c$qcov.fasta
  rm $sequencefile.include.fasta
  cleannumseq=""
  cleannumseq=`grep -v ">" $workdir/$seqname.cleanedOf.$dbname.p$perid.c$qcov.fasta|wc -lc|sed -r 's/\s+/\t/g'|sed -r 's/^\s+//g'`

  ### count the contaminated bases
  perl $scriptfolder/Clean_FastaSeq_of_Numbers-V1.0.pl $sequencefile.removed.fasta ${sequencefile}.removed.fasta.clean
  mv ${sequencefile}.removed.fasta.clean $workdir/$seqname.contaminatedWith.$dbname.p$perid.c$qcov.fasta
  rm $sequencefile.removed.fasta
  contaminationnumseq=""
  contaminationnumseq=`grep -v ">" $workdir/$seqname.contaminatedWith.$dbname.p$perid.c$qcov.fasta|wc -lc|sed -r 's/\s+/\t/g'|sed -r 's/^\s+//g'`

  printf "\n$seqname\t$totalnumseq\t$cleannumseq\t$contaminationnumseq" >> $workdir/$Resultfile

}
#######################################################################################################################################################


while getopts “e:c:s:o:d:p:q:b:f” OPTION
do
     case $OPTION in
         c)
              BacterialGenome=$OPTARG
              ;;
         s)

              sequencefile=$OPTARG
              ;;
         o)
              Resultfile=$OPTARG
              ;;
         d)
              curdir=$OPTARG
              ;;
         p)
              perid=$OPTARG
              ;;
         q)
              qcov=$OPTARG
              ;;
         b)
              noblast=1
              ;;
         f)
              no_format=1
              ;;
         e)
              extra_blast=$OPTARG
              ;;
         ?)
             usage
             exit
             ;;
     esac
done


: ${Resultfile:="Contamination_summary.log"}
: ${curdir:=$(pwd)}
: ${blastdatabase:="/home/Databases/BlastDatabases/Bacteria.Genomes.and.plasmids.ncbi.20130719"}
: ${qcov:=90}
: ${perid:=95}
scriptfolder="/home/ratnesh.singh/Scripts/Perl";



if [ -z $sequencefile ];then
  printf "Usage: Script_name  -options \n
        Options:\n
        -s Sequence file path to clean. [ Required ]\n
        -c contaminating sequence file path [ /home/Databases/BlastDatabases/Bacteria.Genomes.and.plasmids.ncbi.20130719 ]\n
        -p percent identity cutoff to label sequence as contamination [95].
        -q percent of total query length to be matched to call as contamination[90]
        -o output file name to save summary [Contamination_summary.log]\n
        -d save results in dir.[ $(pwd) ]\n
        -b skip blast step.
        -f dont try to create blast database. use the fie name as it is.
        -e extra options passed to blast.place them in quotes.
        "
  exit;
fi




#for file in `ls $curdir|grep "$sequencefile"`; do
#  ### Cal_contamination sequence_file blast_database outputfile
#Cal_contamination "$curdir/$file" $BacterialGenome  "$curdir/$Resultfile.log"
#done


for file in "$sequencefile"; do
  ### Cal_contamination sequence_file blast_database outputfile save_folder
Cal_contamination "$file" "$BacterialGenome"  "$Resultfile.log" "$curdir" "$extra_blast"

done
