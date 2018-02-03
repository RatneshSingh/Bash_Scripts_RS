######tblatx Arabidopsis genes sequences with Cpapaya cds
script_folder="/home/ratnesh.singh/Scripts/Perl"
blast_parser="$script_folder/Blast_parser_BIOPERL_V1.0.pl"
blast_table_sum_redundant="$script_folder/Blast_Table_Sumup_Redundant_hits_V3.0.pl"
search_seq_script="$script_folder/search_seq-V3.4.pl"

if [[ ! -e  $blast_parser  ]]; then
  echo "Unable to find Script:$blast_parser"
  exit
fi
if [[ ! -e  $blast_table_sum_redundant  ]]; then
  echo "Unable to find Script:$blast_table_sum_redundant"
  exit
fi
if [[ ! -e  $search_seq_script  ]]; then
  echo "Unable to find Script:$search_seq_script"
  exit
fi


##############################################################



#############################################################
##############################################################
function usage(){
 echo $1
 echo "USAGE: $0 ...options
 Options: Atleast one query and one database is required
 -s Nucleotide Query Sequence | -S Protein Query Sequence
 -c cds_db
 -p prot_db
 -g genome_db
 -e evalue [1e-10]
 -f final query coverage to define homolog[50]
 -a extra_seq_to_extract_alongwith_hit_region [1000]
 -n output_prefix [Homolog]
 -G Gene_prefix to add in each homologs name
 -b Do not run blast step.
 -B Break and cluster alignments based on tree
 -m prog. Run multiple alignment of homologs with prog
 -d design qRT-PCR primers for tophits.
 -N Dont run blast with NR database for annotation purpose.
 -r Delete old result files.
 "
echo
exit


}


function run_muscle(){
  infile=$1
  outfile=$2
  muscle -in $infile -out $outfile
}

function run_clustalw(){
  infile=$1
  outfile=$2
  clustalw2  -align -infile=$infile -outfile=$outfile -outorder=aligned -QUIET
}

function run_clustalw2(){
  infile=$1
  outfile=$2
  clustalw2  -align -infile=$infile -outfile=$outfile -outorder=aligned -QUIET
}

function run_cobalt(){
  infile=$1
  outfile=$2
  cobalt  -rpsdb /mnt/yulab_16TB/Databases/COBALT_db/cdd_clique_0.75 -i $infile -outfmt mfasta > $outfile
}

function run_blast(){
  local query=$1
  local database=$2
  local outputfile=$3
  local blast_name=$4
  local evalue=$5
  local options=$6

  : ${evalue:="1e-5"}

  if [[ "$blast_name" == "blastn"  ]] ||   [[ "$blast_name" == "tblastx"  ]] ||  [[ "$blast_name" == "tblastn" ]] ;then
    echo ""
    echo " Making nucl blast database from reference file $database"
    echo ""
    [[ ! -s $database.nhr ]]  && makeblastdb -in $database -dbtype nucl

  elif [[ "$blast_name" == "blastp"  ]] ||   [[ "$blast_name" == "blastx" ]] ; then
    echo ""
    echo " Making protein blast database from reference file $database"
    echo ""
    [[ ! -e $database.phr ]] && makeblastdb -in $database -dbtype prot
  fi


  if [[ "$blast_name" == "tblastx"  ]];then
    bopts=" -xdrop_ungap 300 "
  else
    bopts=" -xdrop_gap 300 -xdrop_gap_final 300 "
  fi

  echo "**Running $blast_name on $query against $database"
  echo ""
  #if [[ -e /home/ratnesh.singh/Scripts/Perl/parralel_blast.pl ]]; then
  #  perl /home/ratnesh.singh/Scripts/Perl/parralel_blast.pl -p $blast_name -query $query -db $database -out $outputfile  -a "word_size 9" -outfmt "6 std qcovhsp qcovs qlen slen" -a "culling_limit 2" -a "max_target_seqs 1" -a "num_threads 45" -a "evalue $evalue"  -a "xdrop_gap  300" -a "xdrop_gap_final 300" -n 20
  #else
    cmd="$blast_name -query $query -db $database -out $outputfile -num_threads 20 -evalue $evalue $bopts $options"
    echo "Running command: $cmd"
    $cmd || echo "ERROR:Blast run failed:${cmd}"
  ###-word_size 9 -outfmt '6 std qcovhsp qcovs qlen slen' -culling_limit 2 -max_target_seqs 1 -num_threads 45 -evalue $evalue  -xdrop_gap  300  -xdrop_gap_final 300
  #fi

  echo ""
  echo "Finished running $blast_name"
  echo ""
}


function run_transdecoder(){
  infile=$1
  outfile=$2

  TransDecoder.LongOrfs -t $infile
  blastp -query ${infile}.transdecoder_dir/longest_orfs.pep -db /home/Databases/BlastDatabases/NCBI/refseq_protein out longest_orfs.pep.vs.refseqProt.blastp.table  -out ${infile}.transdecoder_dir/longest_orfs.pep.uniref.blastp.table -num_threads 20 -evalue 1e-1 -max_target_seqs 1 -outfmt 6

  TransDecoder.Predict -t ${infile}  --retain_blastp_hits ${infile}.transdecoder_dir/longest_orfs.pep.uniref.blastp.table --single_best_only



}



function write_groups(){
  grp_file=$1
  seq_file=$2
  perl -e '
  use warnings;
  use strict;
  no autovivification;
      open GRP,$ARGV[0];
      my%list;
      my%sum_grp;
      my%seq;
    print "\nOpening grp file: $ARGV[0]\n";
      $/="\n";
      my $grp=0;
      while(<GRP>){
        s/\s*$|^\s*|^>//g;
        next if m/^\s*$/;
        next if m/^Number\s+of\s+clusters/;
        $grp = $1 if (m/^Cluster\s+([\w\d]+)[\;\s]+size\=\d+|^(unclustered)\s*points/);
        next if (m/^Cluster\s+(\d+)\s\;\s+size\=\d+|^(unclustered)\s*points/);
        $list{$_}=$grp;
        my$sp=$_;
        $sp=~s/_[\S]*//g;
        $sum_grp{$grp}{$sp}+=1;

        #print "\nAdded *$_* in group $grp";
      }
    close(GRP);
    print "\nOpening seq file: $ARGV[1]\n";
    $/="\n>";
      open(SEQ,$ARGV[1]);
      while(<SEQ>){
        my($head,@seq)=split(/\n/,$_);
        $head=~s/\s*$|^\s*|^>//g;
        my$shead=$head;    ## create short headers to match cluster titles.
        $shead=~s/\s+[^\n]*//g;
        #print "\nProcessing seq *$head*";
        #next if ! exists $list{$head};
        my$seq=join("",@seq);
        $seq=~s/[^\w\-]+//g;
        $seq{$list{$shead}}{$head}=$seq;
        #$seq{$list{$head}}{$head}=$seq;

      }
      open OUT2,">",$ARGV[1]."_clust_summary.table";
      for my$group(keys %seq){
        print "\n**********Writing seqs for: $group\n";
        open OUT,">",$ARGV[1]."_clust_".$group.".fasta";
        for my$heads(keys %{$seq{$group}}){
          print OUT ">$heads\n$seq{$group}{$heads}\n";
          #print "\n$heads  for group $group\n";
        }
        close(OUT);
        print OUT2 "\nCluster_ $group";
        for my$sps(sort keys %{$sum_grp{$group}}){
          print OUT2 "\t$sps\:$sum_grp{$group}{$sps}";

        }

      }
      close(OUT2);


    ' $grp_file $seq_file




}



function align_N_tree (){
         sfile=$1
         for aligner in muscle pnpprobs clustalo; do
            if exists ${aligner}; then
                 [[ ! -s ${sfile} ]] || run_${aligner} ${sfile} ${sfile%.*}.${aligner}.fasta
                 [[ ! -s ${sfile%.*}.${aligner}.fasta ]] || make_nj_tree ${sfile%.*}.${aligner}.fasta ${sfile%.*}.${aligner}.fasta.NJ
            fi
            done

}


function get_hmm_homologs(){

  #get_hmm_homologs ${houtname}.prot.fasta_hmmscan.domtblout.pfam.homolog.list $prot_db_fasta prot ${houtname}
          qseq=$1
          domtable=$2
          db=$3
          dbtype=$4
          foutname=$5

          perl $search_seq_script -s  $db  -l ${domtable} -o ${domtable%.*}.${dbtype}.fasta -m exact  -d " " -c 1
          ## outfile name
          hmmout_name=hmmscan_${foutname}.andOthers.${dbtype}
          [[ ! -s ${db} ]] || join_files ${hmmout_name}.fasta ${qseq}  ${domtable%.*}.${dbtype}.fasta

          ## alignments and tree
          for sfile in ${hmmout_name}.fasta ${domtable%.*}.${dbtype}.fasta; do
            align_N_tree $sfile
          done

          ## pnpprobs -> clustalo alignments
          for sfile in ${hmmout_name}.pnpprobs ${domtable%.*}.${dbtype}.pnpprobs; do
            for aligner in clustalo; do
              [[ ! -s ${sfile}.fasta ]] || run_${aligner} ${sfile}.fasta ${sfile}.${aligner}.fasta
              [[ ! -s ${sfile}.${aligner}.fasta ]] || make_nj_tree ${sfile}.${aligner}.fasta ${sfile}.${aligner}.fasta.NJ
            done
          done

        if [[ $break_clust ]];then
            for clfile in ${hmmout_name}.pnpprobs.fasta  ${hmmout_name}.pnpprobs.clustalo.fasta; do
                [[ ! -s $clfile ]] || break_in_clusters $clfile
            done

            for file in $(ls ${hmmout_name}.{pnpprobs,pnpprobs.clustalo}.fasta_clust_*|grep ".fasta$"); do
                sed -i 's/-//g' $file;
                run_pnpprobs ${file} ${file/.fasta/}.pnpprobs.fasta;
                run_clustalo ${file} ${file/.fasta/}.clustalo.fasta;
                run_clustalo ${file/.fasta/}.pnpprobs.fasta ${file/.fasta/}.pnpprobs.clustalo.fasta

                make_nj_tree ${file/.fasta/}.pnpprobs.fasta;
                make_nj_tree ${file/.fasta/}.clustalo.fasta;
                make_nj_tree ${file/.fasta/}.pnpprobs.clustalo.fasta
            done
        fi
          rename ".ph" ".phb" *.ph
}


function run_pnpprobs(){
  infile=$1
  outfile=$2
  if exists pnpprobs; then
    [[ ! -z $outfile ]] || outfile=${infile%.*}.pnpprobs.fasta
    sed -i 's/\*//g' ${infile}
    pnpprobs -num_threads 20 -c 5 -ir 500  --outfile ${outfile} $infile
    run_normd ${outfile}
  else
      echo "pnpprobs is not found. Make sure it is installed and in PATH"

  fi
}

function run_clustalo(){
  infile=$1
  outfile=$2
    if exists clustalo; then
      [[ ! -z $outfile ]] || outfile=${infile%.*}.clustalo.fasta
      sed -i 's/\*//g' ${infile}
      clustalo -i $infile -o $outfile --iter=2 --threads=40 --force
      run_normd ${outfile}
    else
      echo "clustalo is not found. Make sure it is installed and in PATH"

    fi

}

function run_muscle(){
  infile=$1
  outfile=$2
    if exists muscle; then
      [[ ! -z $outfile ]] || outfile=${infile%.*}.muscle.fasta
      sed -i 's/\*//g' ${infile}
      muscle -in $infile -out $outfile -quiet
      run_normd ${outfile}
    else
      echo "muscle is not found. Make sure it is installed and in PATH"

    fi
}


function run_hmmscan(){
  prot_seq=$1
  if exists hmmscan;then
    hmmscan -o ${prot_seq}_hmmscan.out --tblout ${prot_seq}_hmmscan.tblout --pfamtblout ${prot_seq}_hmmscan.pfamtblout --domtblout ${prot_seq}_hmmscan.domtblout --noali --cpu 10 /home/Databases/pfam/Pfam-A.hmm ${prot_seq}
  else
      echo "hmmscan was not found. Make sure it is installed and in PATH"

  fi
}


function make_nj_tree(){
  infile=$1
  outfile=$2
  if exists clustalw2;then
  [[ ! -z $outfile ]] || outfile=${infile}.NJ
  clustalw2 -INFILE=${infile} -TREE -BOOTSTRAP  -QUIET -BOOTLABELS=node -CLUSTERING=NJ -OUTFILE=$outfile || failed=1
  if [[ ! -z ${failed}  && -s ${infile} ]]; then
    perl ~/Scripts/Perl/Clean_FastaSeq_of_Numbers-V1.0.pl ${infile} ${infile%.*}.nname.fasta
    clustalw2 -INFILE=${infile%.*}.nname.fasta -TREE -BOOTSTRAP  -QUIET -BOOTLABELS=node -CLUSTERING=NJ -OUTFILE=$outfile
  fi
    else
      echo "clustalw was not found. Make sure it is installed and in PATH"

  fi
}

function break_in_clusters(){
  infile=$1
  outfile=$2
  [[ ! -z $outfile ]] || outfile=${infile%.*}.cluspack.grp

  [[ ! -s $infile ]] || cluspack ${infile} -dt=alignment -cm=bionj -nbc=secator -oclu=${outfile}
  [[ ! -s $infile ]] || write_groups ${outfile} $infile
}

function cluster_N_realign(){
  infile=$1
  outfile=$2

  [[ ! -z $outfile ]] || outfile=${infile%.*}.cluspack.grp
  break_in_clusters $infile $outfile
  for file in $(ls ${infile}_clust_*|grep ".fasta$"); do
      sed -i 's/-//g' $file;
      run_pnpprobs ${file} ${file/.fasta/}.pnpprobs.fasta;
      run_clustalo ${file} ${file/.fasta/}.clustalo.fasta;
      run_clustalo ${file/.fasta/}.pnpprobs.fasta ${file/.fasta/}.pnpprobs.clustalo.fasta

      make_nj_tree ${file/.fasta/}.pnpprobs.fasta;
      make_nj_tree ${file/.fasta/}.clustalo.fasta;
      make_nj_tree ${file/.fasta/}.pnpprobs.clustalo.fasta
  done

    rename ".ph" ".phb" *.ph



}

function exists()
{
  command -v "$1" >/dev/null 2>&1
}


function run_normd(){
  aln=$1
  if exists normd;then
    normd ${aln} > ${aln}.normd
    sleep 0.01
  else
     echo "normd was not found. Make sure it is installed and in PATH"
  fi
}


function join_files(){
    outfile=$1

    for files in "${@:2}";do
        [[ -e $files ]] && { cat $files; echo; }

    done > $outfile
}

##############################################################





while getopts "s:S:c:p:g:e:a:n:bf:m:t:dNrG:B" OPTION
do
     case $OPTION in
         s)
              query_seq_fasta=$OPTARG
              ;;
         S)
              query_pep_fasta=$OPTARG
              ;;
         c)

              cds_db_fasta=$OPTARG
              ;;
         p)

              prot_db_fasta=$OPTARG
              ;;

         g)
              genome_db_fasta=$OPTARG

              ;;

          e)
              bl_eval=$OPTARG
              ;;
          f)
              fqcm=$OPTARG
              ;;

          a)
              extra_seq=$OPTARG
              ;;

          n)
              file_prefix=$OPTARG
              ;;
          G)
              gene_prefix=$OPTARG
              ;;
          b)
              no_blast=1
              ;;
          B)
              break_clust=1
              ;;
          m)
              multi_align=$OPTARG
              ;;
          t)
              gff=$OPTARG
              ;;
          d)
              design=1
              ;;
          N)
              NoNR=1
              ;;
          r)
              replace=1
              ;;
          ?)
             usage
             exit
             ;;
     esac
done


#: ${cds_db_fasta:="/mnt/ratnesh.singh/Papaya/Papaya_WT_diminutive_RNA_seq_Data/Papaya_All_SMS_DMS_trinity_assembly_NO_DMS-L-0_plus_Cpapaya_113_cds.fasta"}
#: ${prot_db_fasta:="/home/Databases/BlastDatabases/phytozome_data/protein_primaryTranscriptOnly/Cpapaya_113_ASGPBv0.4.protein_primaryTranscriptOnly.fa"}
#: ${genome_db_fasta:="/home/Databases/BlastDatabases/phytozome_data/genomes/Cpapaya_113.fa"}
: ${file_prefix:="Homologs"}
: ${bl_eval:="1e-1"}
## query coverage to be qualified as homolog.
: ${fqcm:=50}

## extract extra sequence from genomic in sequence is prodided.
: ${extra_seq:=1000}

### check for parameters

if [[ -z $query_seq_fasta && -z $query_pep_fasta ]]; then
    usage "No query file detected"
fi


if [[ ! -z $query_seq_fasta && ! -z $query_pep_fasta ]]; then
    usage "Use either cds or protein as query. Not both"
fi


if [[ -z $cds_db_fasta && -z $prot_db_fasta ]]; then
    usage "No database file detected. Please provide atleast one database"
fi



if [[ ! -e $query_seq_fasta ]] && [[ ! -e $query_pep_fasta ]]; then
 echo "Unable to open sequence file: $query_seq_fasta $query_pep_fasta "
 echo
 exit
elif [[ ! -z  $cds_db_fasta && ! -e  $cds_db_fasta ]]; then
 echo "Unable to open cds database file: $cds_db_fasta"
 echo
elif [[ ! -z $prot_db_fasta && ! -e $prot_db_fasta ]]; then
 echo "Unable to open protein database file: $prot_db_fasta"
else
 echo "Searching homolog for $query_seq_fasta $query_pep_fasta in $cds_db_fasta , $prot_db_fasta "
 echo
fi

## define blast type to run

if [[ ! -z $query_seq_fasta ]]; then
    nblast='tblastx'
    pblast='blastx'
    query=$query_seq_fasta

elif [[ ! -z $query_pep_fasta ]]; then
    nblast='tblastn'
    pblast='blastp'
    query=$query_pep_fasta
fi

perl ~/Scripts/Perl/Clean_FastaSeq_of_Numbers-V1.0.pl ${query} ${query%.*}.nname.${query##*.}

query=${query%.*}.nname.${query##*.}

#for file in $query; do
    echo "Looking for homolog of $file in $cds_db_fasta $prot_db_fasta "
    echo "Running blast between $file and $cds_db_fasta $prot_db_fasta ...."

    ### assign variables to complex names for outfiles.
    file=${query##*/}
    pblast_out=${file%.*}_Vs_${prot_db_fasta##*/}.${pblast}
    nblast_out=${file%.*}_Vs_${cds_db_fasta##*/}.${nblast}
    bothblast_out=${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}.blast_result
    houtname=${file_prefix}_of_${file%.*}

    if [[ $replace ]]; then
         rm *${houtname}* *${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}* *${nblast_out}* *${pblast_out}*
    fi


    if [[ ! $no_blast ]]; then
        [[ -e $prot_db_fasta && ! -e $pblast_out ]] && run_blast   $query    $prot_db_fasta   $pblast_out  ${pblast} ${bl_eval}
        [[ -e $cds_db_fasta  && ! -e $nblast_out ]]  && run_blast   $query    $cds_db_fasta    $nblast_out  ${nblast} ${bl_eval}
    fi

    sleep 0.25
    #echo "" > ${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}.blast.temp
    #[[ -s $pblast_out ]] && awk 1  $pblast_out > ${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}.blast.temp
    #[[ -s $nblast_out ]] && awk 1  $nblast_out >> ${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}.blast.temp

    join_files $bothblast_out $pblast_out $nblast_out

    #mv ${file%.*}_Vs_${cds_db_fasta##*/}.N.${prot_db_fasta##*/}.blast.temp $bothblast_out
    #
    echo
    echo "Parsing blast result...."
    perl $blast_parser -i $bothblast_out -o ${bothblast_out}.table -e $bl_eval
    # getting tophit from blast table
    sort -k1,1 -k12,12nr -k11,11n  ${bothblast_out}.table | sort -u -k1,1 --merge > ${bothblast_out}.tophit.table

    echo
    echo "Removing redundant hits from blast result...."
    perl $blast_table_sum_redundant -b ${bothblast_out}.table -a 1000 -o
    #
    echo
    echo "Summarizing blast results...."
    sort -rnk5,5 ${bothblast_out}.table.subject.summary.table|grep -v "Subject"|cut -f1,7>${bothblast_out}.table.subject.summary.table.list
    cut -f 2 ${bothblast_out}.tophit.table | sort|uniq >  ${bothblast_out}.tophit.table.list

    ## getting list of sequences to extract from database.
    sed -i '/^>/s/lcl|//' ${bothblast_out}.table.subject.summary.table.list
    cat ${bothblast_out}.table.subject.summary.table.list|sort|uniq > ${bothblast_out}.table.subject.summary.table.list1
    mv ${bothblast_out}.table.subject.summary.table.list1 ${bothblast_out}.table.subject.summary.table.list

    echo
    echo "Getting sequences homologous to $file from database $cds_db_fasta $prot_db_fasta ...."
    ## for prot
    [[ -s $cds_db_fasta ]]  && perl $search_seq_script -s  $cds_db_fasta  -l ${bothblast_out}.table.subject.summary.table.list -o ${houtname}.cds.fasta -m exact  -d " " -c 1
    [[ -s $cds_db_fasta ]] && perl $search_seq_script -s  $cds_db_fasta  -l ${bothblast_out}.tophit.table.list -o ${houtname}.cds.tophit.fasta -m exact  -d " " -c 1
    ## for cds
    [[ -s $prot_db_fasta  ]] && perl $search_seq_script -s  $prot_db_fasta -l ${bothblast_out}.table.subject.summary.table.list -o ${houtname}.prot.fasta -m exact -d " " -c 1
    [[ -s $prot_db_fasta  ]] && perl $search_seq_script -s  $prot_db_fasta -l ${bothblast_out}.tophit.table.list -o ${houtname}.prot.tophit.fasta -m exact  -d " " -c 1

    #### exit if no names were picked by blast.
    [[ $(cat ${houtname}.{cds,prot}.fasta | grep -c ">" ) < 1 ]] && { echo "No homologs were picked by blast. script will exit now.";exit; }

    ### extract sequences from databases and blast them to NCBI NR for annotation
    if [[ $( grep -c ">" ${houtname}.cds.fasta  ) > 0 ]]; then
      echo
       [[ ! $no_blast ]] && echo "Running blastx for homologous sequences against NCBI protein database locally at evalue 1e-20...."
       [[ ! $NoNR && ! $no_blast ]] && screen -m -d -S NRdb_blastx_run_from_findHomologScript bash -c "blastx -query ${houtname}.cds.fasta -out ${houtname}_VS_refseqdb.blastx.table -db /home/Databases/BlastDatabases/NCBI/refseq_protein -num_threads 20 -evalue 1e-5 -max_target_seqs 5 -culling_limit 2 -outfmt '6 std qcovhsp qcovs qlen slen stitle' -gilist /home/Databases/BlastDatabases/NCBI/viridiplantae_refSeq_gi_list.txt"
      echo ;
      if [[ -s $genome_db_fasta ]]; then
        echo "Running blastn for Papaya homologs against db genome...."
        [[ -s $genome_db_fasta ]] && blastn -query ${houtname}.cds.fasta -db $genome_db_fasta -out ${houtname}.Vs.${genome_db_fasta##*/}.blastn -evalue 1 -num_threads 20
        echo "Parsing Genomic Blastn result...."
        [[ -s ${houtname}.Vs.${genome_db_fasta##*/}.blastn ]]  && perl $blast_parser -i ${houtname}.Vs.${genome_db_fasta##*/}.blastn -o ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table -e 1
        echo "Removing redundant hits and extracting homologous region from papaya genome with $extra_seq extra at each end...."
        [[ -s ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table ]] && perl $blast_table_sum_redundant -b ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table -s $genome_db_fasta -e 1e-1 -l 40 -p 80 -a 100000 -hd $extra_seq -t $extra_seq  -o
        echo " Combining extracted genomic sequences to query sequence to align."
        #[[ -s ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table.seq.fasta ]] && awk 1  ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table.seq.fasta   ${houtname}.cds.fasta >${houtname}.cds.genomic.fasta
        join_files ${houtname}.cds.genomic.fasta   ${houtname}.Vs.${genome_db_fasta##*/}.blastn.table.seq.fasta   ${houtname}.cds.fasta
      fi
    fi

    if [[ $(cat ${houtname}.prot.fasta | grep -c ">" ) > 0 ]]; then
      echo
      [[ ! $no_blast ]] && echo "Running blastp for homologous sequences against NCBI protein database locally at evalue 1e-20...."
      [[ ! $NoNR && ! $no_blast ]] && screen -m -d -S NRdb_blastp_run_from_findHomologScript bash -c "blastp -query ${houtname}.prot.fasta -out ${houtname}_VS_refseqdb.blastp.table -db /home/Databases/BlastDatabases/NCBI/refseq_protein -num_threads 20 -evalue 1e-10 -max_target_seqs 5 -culling_limit 2 -outfmt '6 std qcovhsp qcovs qlen slen stitle' -gilist /home/Databases/BlastDatabases/NCBI/viridiplantae_refSeq_gi_list.txt"

      echo ;
      if [[ -s $genome_db_fasta ]]; then
        echo "Running tblastn for $file homologs against $genome_db_fasta ...."
        [[ -s $genome_db_fasta ]] && tblastn -query ${houtname}.prot.fasta -db $genome_db_fasta -out ${houtname}.Vs.${genome_db_fasta##*/}.tblastn -evalue 1 -num_threads 20
        echo "Parsing Genomic tBlastn result...."
        [[ -s ${houtname}.Vs.${genome_db_fasta##*/}.tblastn ]]  && perl $blast_parser -i ${houtname}.Vs.${genome_db_fasta##*/}.tblastn -o ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table -e 1
        echo "Removing redundant hits and extracting homologous region from papaya genome with $extra_seq extra at each end...."
        [[ -s ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table ]] && perl $blast_table_sum_redundant -b ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table -s $genome_db_fasta -e 1e-5 -l 40 -p 80 -a 100000 -hd $extra_seq -t $extra_seq  -o
        echo " Combining extracted genomic sequences to query sequence to align."
        #[[ -s ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table.seq.fasta ]] && awk 1  ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table.seq.fasta   ${houtname}.prot.fasta >${houtname}.prot.genomic.fasta
        join_files  ${houtname}.prot.genomic.fasta  ${houtname}.Vs.${genome_db_fasta##*/}.tblastn.table.seq.fasta   ${houtname}.prot.fasta
      fi
    fi



    ### Blast sequence to other papaya sequences.
    [[ -z $cds_db_fasta ]] || blastn -query ${houtname}.cds.fasta -db $cds_db_fasta -out ${houtname}.Vs.${cds_db_fasta##*/}.blastn.table -outfmt '6 std qcovhsp qcovs qlen slen'   -xdrop_gap 300 -xdrop_gap_final 300  -num_threads 20 -evalue 0.1
    [[ -z $prot_db_fasta ]] || blastp -query ${houtname}.prot.fasta -db $prot_db_fasta -out ${houtname}.Vs.${prot_db_fasta##*/}.blastn.table -outfmt '6 std qcovhsp qcovs qlen slen'   -xdrop_gap 300 -xdrop_gap_final 300  -num_threads 20 -evalue 0.1

    ### remove "lcl|" from blast hit report. it causes problems in finding self hits
    [[ ! -s ${houtname}.Vs.${cds_db_fasta##*/}.blastn.table ]] || sed -i '/^>/s/lcl|//' ${houtname}.Vs.${cds_db_fasta##*/}.blastn.table
    [[ ! -s ${houtname}.Vs.${prot_db_fasta##*/}.blastn.table ]] || sed -i '/^>/s/lcl|//' ${houtname}.Vs.${prot_db_fasta##*/}.blastn.table
    ### mask non-uniq regions
    [[ ! -s ${houtname}.cds.fasta ]] || perl $script_folder/findUniq-V1.4.pl -b ${houtname}.Vs.${cds_db_fasta##*/}.blastn.table -o ${houtname}.uniqregions.fasta -s ${houtname}.cds.fasta






    ###### do multiple alignment with selected papaya genes and other homologs

    [[ -z $query_seq_fasta ]] || join_files  ${houtname}.andOthers.fasta  ${houtname}.cds.fasta   ${query}
    [[ -z $query_pep_fasta ]] || join_files  ${houtname}.andOthers.fasta  ${houtname}.prot.fasta  ${query}

    [[ -s ${houtname}.andOthers.fasta ]] && align_N_tree ${houtname}.andOthers.fasta
    [[ -s ${houtname}.cds.fasta ]]       && align_N_tree ${houtname}.cds.fasta
    [[ -s ${houtname}.prot.fasta ]]      && align_N_tree ${houtname}.prot.fasta

    if [[ $break_clust ]]; then
        cluster_N_realign  ${houtname}.andOthers.pnpprobs.fasta;
        cluster_N_realign  ${houtname}.andOthers.clustalo.fasta;
        cluster_N_realign  ${houtname}.andOthers.pnpprobs.clustalo.fasta
    fi


    ### design primer for pRT-PCR for top hits.
    if [[ $design ]];then

      python ~/Scripts/Python/primer3_online_design_mechanize.py  ${houtname}.cds.tophit.fasta

    fi



  if [[ ! -z $query ]]; then
      echo "Running pfam on protein query to identify the domains"
      [[ -s ${query}_hmmscan.tblout ]] || run_hmmscan ${query}
      [[ -s ${houtname}.prot.fasta_hmmscan.tblout ]] || run_hmmscan ${houtname}.prot.fasta
  fi

 ## Parse pfam domtblout to identify homologous sequences based on sahred domains.
  if [[ -s ${houtname}.prot.fasta_hmmscan.domtblout ]]; then
          python3  ~/Scripts/Python/find_similar_pfamdomain.py --query ${query}_hmmscan.domtblout --homologs ${houtname}.prot.fasta_hmmscan.domtblout --out ${houtname}.prot.fasta_hmmscan.domtblout.pfam.homolog.list --perc_dom 70 --ival 1e-5
          [[ -z $prot_db_fasta ]] || get_hmm_homologs ${query} ${houtname}.prot.fasta_hmmscan.domtblout.pfam.homolog.list $prot_db_fasta prot ${houtname}
          [[ -z $cds_db_fasta ]]  || get_hmm_homologs $empty ${houtname}.cds.fasta_hmmscan.domtblout.pfam.homolog.list $cds_db_fasta cds ${houtname}
  fi

  if [[ $gene_prefix ]]; then
        for rnfile in ${houtname}.prot.fasta_hmmscan.domtblout.pfam.homolog.cds.fasta ${houtname}.prot.fasta_hmmscan.domtblout.pfam.homolog.prot.fasta ${houtname}.cds.fasta ${houtname}.prot.fasta ${houtname}.cds.tophit.fasta ${houtname}.prot.tophit.fasta; do
            [[ -s $rnfile ]] && sed -i "/^>/s/>/>${gene_prefix}/" $rnfile
        done
  fi




for rmd in *.normd;do
  normd=$(cat $rmd)
  printf "\n${rmd}\t${normd}"
done > Alignments_normd_summary.table



    echo
    echo
#done
