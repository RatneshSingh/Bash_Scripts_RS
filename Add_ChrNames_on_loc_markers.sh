#!/bin/bash

###### parsing options
######## usage function
function usage()
{
cat << EOF
usage: $0 options
Required options:
    -j Joinmap's .loc file
optional
    -w       : working directory where you want all the output from this script should be saved [.]
    -r       : Full path length of Reference Genome to be used for chromosome mapping estimation
    -c       : Catalog.tags.tsv file produced by Stacks
    -m       : markers types to be seperated in seperate files.
             : e.g. nnxnp, lmxll, "nnxnp lmxll abxac" etc.["nnxnp lmxll"]
    -h       : Print usage help
";
EOF
}

##############################################################################################################################
function run_blast(){
  local query=$1
  local database=$2
  local outputfile=$3
  local blast_name=$4

  if [ ! -e $database.nhr ];then
    echo " Making blast database from reference file $database"
    makeblastdb -in $database -dbtype nucl
  fi

  echo "Running $blast_name on $query against $database"
  $blast_name -query $query -db $database -out $outputfile  -word_size 11 -outfmt 6 -max_target_seqs 1 -num_threads 30 -evalue 0.00000001  -xdrop_gap  300  -xdrop_gap_final 300
  echo "Finished running $blast_name"
}







#################################################################################################################################
###### Main program
#################################################################################################################################
while getopts “w:j:r:c:h:m:” OPTION
do
     case $OPTION in
        w)
          wordir=$OPTARG
          ;;

         j)
              jmpfile=$OPTARG
              ;;
         r)

              reference_genome=$OPTARG
              ;;
         c)
              catalog_file=$OPTARG
              ;;
         h)
             usage
             exit
             ;;
        m)
            markers_types=$OPTARG
            ;;

        ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $jmpfile ]]
then
     printf "\n\nERROR: Joimmap .loc files are required to run the script.\n\n"
     usage
     exit 1
fi

####==================================================================================================
######################################################################################
###### extract loci reads from joinmap file and create nnxnp and lmxll files.
###function process_loc_file(){

curdir=$(pwd)
jmp=$(realpath $jmpfile)
jmpfile=$jmp
jmpdir=${jmp%/*}
#wordir=${jmpdir%/*}
  : ${markers_types:="nnxnp lmxll"}
  : ${catalog_file:=batch_1.catalog.tags.tsv}


if [[ -e $reference_genome ]]; then
      reference_file=${reference_genome##*/}
      echo "Working for Joinmapfile:$jmpfile"
      echo "$jmpfile: Get marker sequences from $catalog_file"
      cut -f1 $jmpfile|sed  '1,/nind/d'|awk '/individual names:/{seen=1} seen {exit}; {print}'|sort|grep -v "^$" >$jmpfile.catlog.list
      cut -f3,10 $catalog_file|sort -k1n,1 > $catalog_file.catalog_seqs.table
      echo "$jmpfile: Create fasta file to be used as query for blast"

      ##### Takes too long to create fasta. perl hash (below) is much faster for this purpose.
      perl -e '
                open INFILE,"$ARGV[0]";
                while(<INFILE>){
                  $_=~s/^\s+//g;
                  ($header,$sequence)=split /\s+/,$_;
                  $header=~s/^\s+|\s+$//g;
                  $sequence=~s/^\s+|\s+$//g;
                  $seq_hash{$header}=$sequence;
                }
      open LIST,"$ARGV[1]";
      open FASTA,">$ARGV[2]";
      while(<LIST>){

        $_=~s/^\s+|\s+$//g;
        print FASTA ">$_\n$seq_hash{$_}\n";

      }
      close FASTA;
      close LIST;
      close INFILE;
      ' $catalog_file.catalog_seqs.table $jmpfile.catlog.list $jmpfile.catlog.fasta



      echo "$jmpfile: Step: Run blastn to get chromosome coordinates for marker"
      run_blast   $jmpfile.catlog.fasta  $reference_genome  $jmpfile.catlog.vs.$reference_file.blastn.table  blastn

      ## sort and merge multiple hits for each query. select best hit only
      echo "$jmpfile:  Step: Parse blastn result to get best hit"
      sort -k1,1 -k12,12nr -k11,11n  $jmpfile.catlog.vs.$reference_file.blastn.table | sort -u -k1,1 --merge|cut -f1,2,9|awk '{print $1"\t"$1"_"$2"_"$3}' >$jmpfile.catlog.vs.$reference_file.blastn.besthit_conversion.table
      echo "$jmpfile:  Step: Adding coordinates to marker names "
      awk 'NR==1 { next } FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' $jmpfile.catlog.vs.$reference_file.blastn.besthit_conversion.table $jmpfile >$jmpfile.wrt_${reference_file:0:5}.loc
      sed -i 's/scaffold_/S/' $jmpfile.wrt_${reference_file:0:5}.loc
      sed -i 's/Qyu_RAD_index/QRi/' $jmpfile.wrt_${reference_file:0:5}.loc
      sed -i 's/L003_R1_001/P1/' $jmpfile.wrt_${reference_file:0:5}.loc
      ##### select name till first "|" sign.
      sed -i 's/|[^ tab]\+//' $jmpfile.wrt_${reference_file:0:5}.loc
fi


    locfile=$jmpfile.wrt_${reference_file:0:5}.loc

    [[ ! -e $reference_genome ]] && locfile=$jmpfile;

    ###### filter joinmap file to get nnxnp and lmxll loci out in seperate files.
    #[[ $(cd $stack_output_folder) ]] || exit_on_error "cd $stack_output_folder"
    #for file in `ls|grep -e "loc$"|egrep -v "(nnxnp|lmxll)|grep modname.loc"`; do
    #cd $wordir
     for file in ${locfile##*/} ; do
      for gentype in $markers_types; do

        echo "Step: Seperating $gentype from Joinmap file and saving as new file: $gentype.$file "
          nloc=$(grep -c $gentype $locfile)
          name=${gentype}_$(grep "^name =" $locfile|awk '{print $3}')
          #name=$(grep "^name =" $file|awk '{print $3}')
          #name=${gentype}_$name
          name=${name:0:20}
          printf "name = $name\n">$gentype.${file##*/}
          egrep "(^popt =)" $locfile>> $gentype.${file##*/}
          printf "nloc = $nloc\n" >>$gentype.${file##*/}
          egrep "(^nind =)" $locfile>> $gentype.${file##*/}
          printf "\n" >>$gentype.${file##*/}
          echo >>$gentype.${file##*/}
          grep "<$gentype>" $locfile>>$gentype.${file##*/}
          echo >>$gentype.${file##*/}
          awk '/individual names:/ {seen = 1}; seen {print}' $locfile >>$gentype.${file##*/}
      done
    done
######}

cd $curdir;
exit