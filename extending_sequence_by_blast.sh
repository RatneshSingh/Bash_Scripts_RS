########Variables
read_database=$1
seq2extend=$2
cycle=$3
current_folder=`pwd`
########
${cycle:=3}
db_name=`echo ${read_database##*/}`
db_name=`echo ${db_name%%.*}`
seqname=`echo ${seq2extend%%.*}`
############################################################
function extend_seq(){
  local seq2extend=$1
  local read_database=$2
  local round=$3

  local db_name=`echo ${read_database##*/}`
  local db_name=`echo ${db_name%%.*}`
  local seqname=`echo ${seq2extend%%.*}`
  local count=0;

  ln -s  $seq2extend $seqname.Vs.$db_name.round0.blastn.table.list.fasta


  for((c=1;c<=$round;c++)); do
    let d=c-1;
    blastn -query $seqname.Vs.$db_name.round$d.blastn.table.list.fasta -db $read_database -out $seqname.Vs.$db_name.round$c.blastn -num_threads 30
    perl ~/Scripts/Perl/Blast_parser_BIOPERL_V1.0.pl -i $seqname.Vs.$db_name.round$c.blastn -o $seqname.Vs.$db_name.round$c.blastn.table -p 95
    cut -f2 $seqname.Vs.$db_name.round$c.blastn.table|sort|uniq >$seqname.Vs.$db_name.round$c.blastn.table.list
    perl ~/Scripts/Perl/search_seq-V3.0.pl -s $read_database -l $seqname.Vs.$db_name.round$c.blastn.table.list -o $seqname.Vs.$db_name.round$c.blastn.table.list.fasta -d " " -c 1
    rm $seqname.final_round.fasta
    ln -s $seqname.Vs.$db_name.round$c.blastn.table.list.fasta $seqname.final_round.fasta

  done

  perl -ne '$/="\n>"; ($head,@seq)=split /\n/,$_; $seq=join("",@seq);$seq=~s/>//g; $head=~s/>//;$count=0;foreach $seq_s(split(/N+/,$seq)){print "\n>$head\_$count\n$seq_s";$count++}' $seqname.final_round.fasta >$seqname.final_round.fasta.fragmentedAtN.fasta
  cap3 $seqname.final_round.fasta.fragmentedAtN.fasta


}
###########################################################
cd $current_folder/
extend_seq $seq2extend $read_database $cycle
