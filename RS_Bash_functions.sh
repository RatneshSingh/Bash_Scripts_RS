function CleanBactCont(){
    #### usage CleanBactCont  1:sequencefile   2:percent_id[95]   3:query_coverege[80]   4:alignment_length[100]    5:merge_hits/no_merge[merge_hits]   6:allowed_length_to_merge_hits[100]
    local sequencefile=$1
    local        perid=$2
    local         qcov=$3
    local          len=$4
    local        merge=$5
    local      allowed=$6
    scriptfolder="/home/ratnesh.singh/Scripts/Perl";
    BacterialGenome="/home/Databases/BlastDatabases/Bacteria.plasmids.ncbi.20130719.fna"

    ####Set the default for variables
    : ${perid:=95}
    : ${qcov:=80}
    : ${len:=100}
    : ${merge:="merge_hits"}
    : ${allowed:=100}


    #### Run blast against Bacterial genome database
    echo;echo;echo
    echo "Running blast againt Bacterial genome database"
    blastn -query $sequencefile -db $BacterialGenome -out $sequencefile.BacterialGenome.blastn -num_threads 20
    echo;echo;echo;echo "Converting blast result into table format"
    perl $scriptfolder/Blast_parser_BIOPERL_V1.0.pl -i $sequencefile.BacterialGenome.blastn -o $sequencefile.BacterialGenome.blastn.table


  if [ $merge == "no_merge" ]; then
    #### filter directly from blast table. Multiple hits are not merged.
    echo;echo;echo;echo "filter directly from blast table. Multiple hits are not merged."
    cut -f1 $sequencefile.BacterialGenome.blastn.table|sort|uniq > $sequencefile.BacterialGenome.blastn.p$perid.c$qcov.list

  else
    ####filter after merging hits.
    echo;echo;echo;echo "Filtering Blast Result "
    perl $scriptfolder/Blast_Table_Sumup_Redundant_hits_V1.4.pl -b $sequencefile.BacterialGenome.blastn.table -k query -a $allowed -l $len -p $perid -fqcm $qcov|egrep -v "(Printing|Reading|Done|Subject|Total)"|cut -f 1|sort|uniq>$sequencefile.BacterialGenome.blastn.p$perid.c$qcov.list
  fi



    #### Removing contamination from assembled sequence
    echo;echo;echo;echo "Removing contamination from assembled sequence"
    perl $scriptfolder/Remove_selectedSequences-V1.3.pl -s $sequencefile -l $sequencefile.BacterialGenome.blastn.p$perid.c$qcov.list -m match

    #### Reblast filtered sequences to see if any contamination is left.
    echo;echo;echo;echo "Running blast on filtered sequences to see if any contamination is left"
    blastn -query $sequencefile.include.fasta -db $BacterialGenome -out $sequencefile.include.fasta.BacterialGenome.blastn -num_threads 20
    perl $scriptfolder/Blast_parser_BIOPERL_V1.0.pl -i $sequencefile.include.fasta.BacterialGenome.blastn -o $sequencefile.include.fasta.BacterialGenome.blastn.table


    #### blast removed sequences to see what is contamination type
    echo;echo;echo;echo "Running blast on removed sequences to see what is contamination type"
    blastn -query $sequencefile.removed.fasta -db $BacterialGenome -out $sequencefile.removed.fasta.BacterialGenome.blastn -num_threads 20
    perl $scriptfolder/Blast_parser_BIOPERL_V1.0.pl -i $sequencefile.removed.fasta.BacterialGenome.blastn -o $sequencefile.removed.fasta.BacterialGenome.blastn.table
}

#### Remove vector contamination from sequence ends.
#### Usage remove_vector sequencefile vectorfile
#### outputs
function remove_vector(){
    local sequencefile=$1
    local   vectorfile=$2
    local        perid=$3
    local     from_end=$4

    ####Set the default for variables
    : ${perid:=95}
    : ${from_end=1500}


    if [ ! -s $vectorfile.nhr ]; then
      makeblastdb -in $vectorfile -dbtype nucl
    fi

    blastn -query $sequencefile -db $vectorfile -out $sequencefile.Vector.blastn -num_threads 20
    perl $scriptfolder/Blast_parser_BIOPERL_V1.0.pl -i $sequencefile.Vector.blastn -o $sequencefile.Vector.blastn.table -p $perid -b query -m $from_end
    perl $scriptfolder/Blast_Table_Sumup_Redundant_hits_V1.3.pl -b $sequencefile.Vector.blastn.table -a 100 -vectrim -fe 500 -s $sequencefile -k query
}


function Cal_contamination(){

  ####usage
  ####   Cal_contamination 1:sequencefile 2:BacterialGenome 3:Resultfile 4:perid 5:qcov 6:Aln_len 7:merge_hit/no_merge 8:allowed_len_for_merge
  local    sequencefile=$1
  local BacterialGenome=$2
  local      Resultfile=$3
  local           perid=$4
  local            qcov=$5
  local             len=$6
  local           merge=$7
  local         allowed=$8
  local    scriptfolder="/home/ratnesh.singh/Scripts/Perl";


    ####Set the default for variables
    : ${perid:=95}
    : ${qcov:=80}
    : ${len:=100}
    : ${merge:="merge_hits"}
    : ${allowed:=100}
    : ${BacterialGenome:="/home/Databases/BlastDatabases/Bacteria.plasmids.ncbi.20130719.fna"}
  #### usage CleanBactCont  1:sequencefile   2:percent_id[95]   3:query_coverege[80]   4:alignment_length[100]    5:merge_hits/no_merge[merge_hits]   6:allowed_length_to_merge_hits[100]
  CleanBactCont $sequencefile $perid $qcov $len "merge_hits" $allowed
  totalnumseq=""
  totalnumseq=`grep -v ">" $sequencefile.include.fasta|sed -r s/\s+//g|wc`

  contaminationnumseq=""
  contaminationnumseq=`grep -v ">" $sequencefile.removed.fasta|sed -r s/\s+//g|wc`

  tempseqname=`echo $sequencefile|rev|cut -d "/" -f 1|rev`

  echo $tempseqname $totalnumseq $contaminationnumseq >> $Resultfile
}

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
  $blast_name -query $query -db $database -out $outputfile  -word_size 11 -outfmt 6 -max_target_seqs 1 -num_threads 30 -evalue 0.000001  -xdrop_gap  300  -xdrop_gap_final 300
  echo "Finished running $blast_name"
}






function export_stack_data(){
  #### Usage export_stack_data   1:database_name    2:output_folder     3:output_prefix     4:haplo/depth[depth]    5:batch[1]     6:snps_l[1]     7:mark['aa/bb']    8:alle_l[1]
  local database_name=$1
  local output_folder=$2
  local output_prefix=$3
  local haplo=$4
  local batch=$5
  local snps_l=$6
  local mark=$7
  local alle_l=$8
  local hapgeno=$9

  : ${batch:=1}
  : ${output_prefix:="$database_name.export"}
  : ${output_folder:=`pwd`}
  : ${snps_l:=1}
  : ${mark:='aa/bb'}
  : ${alle_l:=1}
  : ${haplo:='depth'}
  : ${hapgeno:='haplo'}

##### Export result from mysql to file
cd $output_folder
for outfmt in "tsv" "xls"; do
if [ $haplo == "haplo" ]; then
  export_sql.pl -D $database_name -b $batch -a $hapgeno -f $output_folder/$output_prefix.$haplo.$outfmt -o $outfmt -F snps_l=1 -F mark=$mark1/$mark2 -F alle_l=1
else
  export_sql.pl -D $database_name -b $batch -a $hapgeno -f $output_folder/$output_prefix.$haplo.$outfmt -o $outfmt -d  -F snps_l=1 -F mark=$mark1/$mark2 -F alle_l=1
fi
done
}





#### Filter stacks result based on coverage and match depth data with the fileterd result.
function filter_stacks(){
  local haplo_stack_file=$1;
  local depth_stack_file=$2;
  local present_in_samples=$3
  : ${present_in_samples:=80}

  ####create temp env variable to send value of filter inside perl script.
  export temp_percent=$present_in_samples
  selected_lines=`
  perl -ne '
    chomp;
    $firstline=$_ if $line1 ne "true";
    $line1="true";
    next if $_=~/^\#/;
    my@header=split(/\t/,$firstline);
    my@data=split(/\t/,$_);

    s/^\s+|\s+$//g foreach(@header);
    s/^\s+|\s+$//g foreach(@data);

    for my $i (0 .. $#header) {$data{$header[$i]}{$data[0]}=$data[$i];}
    $empty=0;
    for my $i (12 .. $#header){
      $data[$i]=~s/\s+//g;
      #print "\t".$header[$i].":".$data[$i];
      #$empty++ if $data[$i] eq "";

      # to exclude index 14 from calculations
      $empty++ if $data[$i] eq "" && $header[$i] ne 'index14';

    }
    $per_cov=eval(($#header-11-$empty)*100/($#header-11));
    next if $per_cov<$ENV{'temp_percent'};
    #TO INCLUDE CARRIZO AS JAPONICA
    $japonica=join("",$data{index5}{$data[0]},$data{index12}{$data[0]},$data{index15}{$data[0]},$data{index18}{$data[0]});
    #TO EXCLUDE INDEX5(CARRIZO) AS JAPONICA
    #$japonica=join("",$data{index12}{$data[0]},$data{index15}{$data[0]},$data{index18}{$data[0]});
    $matrella=join("",$data{index6}{$data[0]},$data{index13}{$data[0]},$data{index19}{$data[0]});
    ####print "\n*******$data[0]:japonica:$japonica\tmatrella:$matrella";
    #### clean, sort and compare snp string to check if any base is common in both.
    $matrella=~s/\W+//g;
    $matrella=join ("", sort (split(//, $matrella)));
    $matrella=~s/(.*)\1/$1/g;

    $japonica=~s/\W+//g;
    $japonica=join ("", sort (split(//,$japonica)));
    $japonica=~s/(.*)\1/$1/g;
    ####print "\t\t\tjaponica:$japonica\tmatrella:$matrella";
    $match='false';
    foreach $jap(split //,$matrella){
      foreach $mat(split //,$japonica){
         $match='true' if lc$jap eq lc$mat;
      }

    }

    ## print the data only if no snp is shared.
    print "\n$data{$header[$0]}{$data[0]}" if  $match eq 'false';
    ;


  ' < $haplo_stack_file
  `
  #### unset env variable.
  unset temp_percent
  #echo $selected_lines
  grep Catalog $haplo_stack_file|head -n 1
  for catalog in $selected_lines; do
    #echo "For catalog:$catalog";
    grep -P "^$catalog\s" $haplo_stack_file;
    grep -P "^$catalog\s" $depth_stack_file;
    #echo $catalog;
  done


}

function assemble_stacks_result(){

   #perl -ne '@array=split /\s+/; $array[6]=~s/\D//g;$base=substr($array[2],$array[6],1); substr($array[2],$array[6],1)=lc$base; print ">\Sequence_".$i++."\n$array[2]\n";' $1>$1.mod.fasta
   cap3 $1.mod.fasta

}

  ######Convert and sort sam file to BAM format.
function sam_to_bam(){
   #### usage: sam_to_bam  sam_file   output_file
  local sam_file=$1
  local output_file=$2
  : ${output_file:=`echo $sam_file|sed 's/.sam/.sorted.bam/'`}

  if [ ! -s $output_file ]; then
  echo "Samtools: conveting and sorting sam alignment to bam format"
    nice samtools view -bS  $sam_file|samtools sort - $output_file
  fi
}


function bam_to_vcf(){
  #### usage: bam_to_vcf  reference_file   bam_file   output_file
  local reference_file=$1      ## reference sequence file to build index from
  local bam_file=$2   ## sequenc file to align against reference folder
  local output_file=$3   ##save output from all the files
  : ${output_file:=`echo $bam_file|sed 's/.bam/.vcf/'`}


  if [ ! -s $output_file.vcf ]; then
    echo "Samtools: Creating pileup file from ${seq_name}.alignedon_${ref_name}.bowtie2.sorted.bam"
    nice samtools mpileup -uf $reference_file $bam_file | bcftools view -vcg - > $output_file
  else
    echo "$output_file exists. Remove or rename  to run mpileup again. "
  fi

}



function run_bowtie(){
  #### usage: run_bowtie  reference_file   sequence_file   output_folder
  local reference_file=$1      ## reference sequence file to build index from
  local sequence_file=$2   ## sequenc file to align against reference folder
  local output_folder=$3   ##save output from all the files
  local sam_to_bam=$4
  local bam_to_bcf=$5
        : ${output_folder:="."}
  ### get the file name from reference file and sequence files in case user provided full path
  local ref_array=(${reference_file//\// })
  local ref_name=${ref_array[${#ref_array[@]}-1]}

  local seq_array=(${sequence_file//\// })
  local seq_name=${seq_array[${#seq_array[@]}-1]}


  mkdir -p $output_folder
  cd       $output_folder
  if [ ! -s "$reference_file.1.ebwt" ];then
    echo "Creating bowtie index file"

    nice bowtie-build $reference_file $reference_file
  else
    echo "bowtie index file for $reference_file already exists in $current_folder. Please remove or rename it to create new index."
  fi

  if [ ! -s $reference_file.fai ]; then
    echo "Creating fai index file"
    perl -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s+|>//g;$sequence =~ s/(.{0,80})/$1\n/g;print ">$header\n$sequence\n";} ' <$reference_file>$reference_file.mod
    mv $reference_file.mod $reference_file
    nice samtools faidx $reference_file
   else
    echo "fai index file for $reference_file already exists in $current_folder. Please remove or rename it to create new index."
  fi

  ##### align Reads on the reference genome using bowtie for SNP identification
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bowtie.sam ];then
    echo "Bowtie: Aligning Reads from $sequence_file on $reference_file"
    nice bowtie --sam -p 40 -n 3 -l 18 --maxbts 800 -y --best --chunkmbs 10000  $reference_file $sequence_file ${seq_name}.alignedon_${ref_name}.bowtie.sam  &> $sequence_file ${seq_name}.alignedon_${ref_name}.bowtie.log
  else
    echo "sam file with the name ${seq_name}.alignedon_${ref_name}.bowtie.sam exists. please rename or delete that file before running bowtie"
  fi
  ######Convert and sort sam file to BAM format.
if [ $sam_to_bam ];then
    if [ ! -s ${seq_name}.alignedon_${ref_name}.bowtie.sorted.bam ]; then
    echo "Samtools: conveting and sorting sam alignment to bam format"
    nice samtools view -bS  ${seq_name}.alignedon_${ref_name}.bowtie.sam|samtools sort - ${seq_name}.alignedon_${ref_name}.bowtie.sorted
  fi
fi

if [ $bam_to_vcf ];then
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bowtie.sorted.pileup ]; then
    echo "Samtools: Creating pileup file from ${seq_name}.alignedon_${ref_name}.bowtie.sorted.bam"
    nice samtools mpileup -uf $reference_file ${seq_name}.alignedon_${ref_name}.bowtie.sorted.bam | bcftools view -vcg - > ${seq_name}.alignedon_${ref_name}.bowtie.sorted.vcf
  else
    echo "${seq_name}.alignedon_${ref_name}.bowtie.sorted.pileup exists. Remove or rename  to run pileup again. "
  fi
fi
}


function bam_to_SNP(){
  local ref_seq=$1
  local out_put_prefix=$2
  local depth=$3
  local sorted_bam_list=$4
  : ${depth:=10}
  ###### find variant in sorted BAM file using SAMtools/bcf tools
  if [ ! -s $out_put_prefix.raw.bcf.varfiltD${depth}.vcf ];then
    echo "Samtools: running mpile on sorted bam alignment."
    nice samtools mpileup -uf $ref_seq $sorted_bam_list | bcftools view -bvcg - > $out_put_prefix.raw.bcf
    nice bcftools view $out_put_prefix.raw.bcf | vcfutils.pl varFilter -D${depth} > $out_put_prefix.raw.bcf.varfiltD${depth}.vcf
  else
  printf "VCF file already exists. "

  fi

}


function run_bowtie2(){
  #### usage: run_bowtie  reference_file   sequence_file   output_folder
  local reference_file=$1      ## reference sequence file to build index from
  local sequence_file=$2   ## sequenc file to align against reference folder
  local output_folder=$3   ##save output from all the files
  local sam_to_bam=$4
  local bam_to_bcf=$5
        : ${output_folder:="."}
        #: ${sam_to_bam:="no"}
        #: ${bam_to_bcf:="no"}
  ### get the file name from reference file and sequence files in case user provided full path
  local ref_array=(${reference_file//\// })
  local ref_name=${ref_array[${#ref_array[@]}-1]}

  local seq_array=(${sequence_file//\// })
  local seq_name=${seq_array[${#seq_array[@]}-1]}


  mkdir -p $output_folder
  cd       $output_folder
  if [ ! -s "$reference_file.1.bt2" ];then
    echo "Creating bowtie2 index file"

    nice bowtie2-build $reference_file $reference_file
  else
    echo "bowtie2 index file for $reference_file already exists in $current_folder. Please remove or rename it to create new index."
  fi

  if [ ! -s $reference_file.fai ]; then
    echo "Creating fai index file"
    perl -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s+|>//g;$sequence =~ s/(.{0,80})/$1\n/g;print ">$header\n$sequence\n";} ' <$reference_file>$reference_file.mod
    mv $reference_file.mod $reference_file
    nice samtools faidx $reference_file
   else
    echo "fai index file for $reference_file already exists in $current_folder. Please remove or rename it to create new index."
  fi

  ##### align Reads on the reference genome using bowtie2 for SNP identification
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bowtie2.sam ];then
    echo "Bowtie2: Aligning Reads from $sequence_file on $reference_file"
    nice bowtie2  -5 20 --phred33  -L 15 -N 1 -p 20 -x $reference_file -U $sequence_file -S ${seq_name}.alignedon_${ref_name}.bowtie2.sam  &> $sequence_file ${seq_name}.alignedon_${ref_name}.bowtie2.log
  else
    echo "sam file with the name ${seq_name}.alignedon_${ref_name}.bowtie2.sam exists. please rename or delete that file before running bowtie2"
  fi


}






function run_novoalign(){
  #### usage: run_bowtie  reference_file   sequence_file   output_folder
  local reference_file=$1      ## reference sequence file to build index from
  local sequence_file=$2   ## sequenc file to align against reference folder
  local output_folder=$3   ##save output from all the files
  local sam_to_bam=$4

        : ${output_folder:="."}
        #: ${sam_to_bam:="no"}
        #: ${bam_to_bcf:="no"}
  ### get the file name from reference file and sequence files in case user provided full path
  local ref_array=(${reference_file//\// })
  local ref_name=${ref_array[${#ref_array[@]}-1]}

  local seq_array=(${sequence_file//\// })
  local seq_name=${seq_array[${#seq_array[@]}-1]}


  mkdir -p $output_folder
  cd       $output_folder
  if [ ! -s "$reference_file.ndx" ];then
    echo "Creating novoindex for $reference_file"

    nice novoindex -k 14 -s 1 -t 30 $reference_file.ndx $reference_file
  else
    echo "novo index file for $reference_file already exists in $current_folder. Please remove or rename it to create new index."
  fi

  ##### align Reads on the reference genome using bowtie2 for SNP identification
  if [ ! -s ${seq_name}.alignedon_${ref_name}.novoalign.sam ];then
    echo "novoalign: Aligning Reads from $sequence_file on $reference_file"
    nice novoalign  -H -x 3 -a -o SAM  -d $reference_file.ndx -f $sequence_file  > ${seq_name}.alignedon_${ref_name}.novoalign.sam  2> ${seq_name}.alignedon_${ref_name}.novoalign.log
  else
    echo "sam file with the name ${seq_name}.alignedon_${ref_name}.novoalign.sam exists. please rename or delete that file before running novoalign"
  fi

  if [ "$sam_to_bam" == "bam" ]; then
    sam_to_bam ${seq_name}.alignedon_${ref_name}.novoalign.sam  $reference_file
  elif [ "$sam_to_bam" == "vcf" ]; then
    sam_to_bam ${seq_name}.alignedon_${ref_name}.novoalign.sam  $reference_file  "vcf"
  fi


}











function run_bwa_samse(){
  #### usage: run_bowtie  reference_file   sequence_file   output_folder
  local reference_file=$1      ## reference sequence file to build index from
  local sequence_file=$2   ## sequenc file to align against reference folder
  local output_folder=$3   ##save output from all the files
  local sam_to_bam=$4
  local bam_to_bcf=$5
        : ${output_folder:="."}
        #: ${sam_to_bam:="no"}
        #: ${bam_to_bcf:="no"}
  ### get the file name from reference file and sequence files in case user provided full path
  local ref_array=(${reference_file//\// })
  local ref_name=${ref_array[${#ref_array[@]}-1]}

  local seq_array=(${sequence_file//\// })
  local seq_name=${seq_array[${#seq_array[@]}-1]}


  mkdir -p $output_folder
  cd       $output_folder
  if [ ! -s "$reference_file.sa" ];then
    echo "Creating bwa index file"

    nice bwa index $reference_file
  else
    echo "bwa index file for $reference_file already exists. Please remove or rename it to create new index."
  fi


  ##### align Reads on the reference genome using bwa for SNP identification
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bwa.sam ];then
    echo "bwa: Aligning Reads from $sequence_file on $reference_file"
    nice bwa aln -l 20 -t 20 $reference_file $sequence_file > ${seq_name}.alignedon_${ref_name}.bwa.aln
    nice bwa samse -L 15 -N 1 -p 20 $reference_file ${seq_name}.alignedon_${ref_name}.bwa.aln $sequence_file >${seq_name}.alignedon_${ref_name}.bwa.sam
  else
    echo "sam file with the name ${seq_name}.alignedon_${ref_name}.bwa.sam exists. please rename or delete that file before running bwa"
  fi

  ######Convert and sort sam file to BAM format.
if [ $sam_to_bam ];then
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bwa.sorted.bam ]; then
  echo "Samtools: conveting and sorting sam alignment to bam format"
    nice samtools view -bS  ${seq_name}.alignedon_${ref_name}.bwa.sam|samtools sort - ${seq_name}.alignedon_${ref_name}.bwa.sorted
  fi
fi

if [ $bam_to_vcf ];then
  if [ ! -s ${seq_name}.alignedon_${ref_name}.bwa.sorted.pileup ]; then
    echo "Samtools: Creating pileup file from ${seq_name}.alignedon_${ref_name}.bwa.sorted.bam"
    nice samtools mpileup -uf $reference_file ${seq_name}.alignedon_${ref_name}.bwa.sorted.bam | bcftools view -bvcg - > ${seq_name}.alignedon_${ref_name}.bwa.sorted.pileup
  else
    echo "${seq_name}.alignedon_${ref_name}.bwa.sorted.pileup exists. Remove or rename  to run pileup again. "
  fi
fi

}

function sam_to_bam(){

  local samfile=$1
  local reference_file=$2
  local bam_to_vcf=$3  #### vcf

  ######Convert and sort sam file to BAM format.

    if [ ! -s ${samfile%sam}sorted.bam ]; then
    echo "Samtools: conveting and sorting sam alignment to bam format"
    nice samtools view -bS  $samfile|samtools sort - ${samfile%sam}sorted
  fi


if [ $bam_to_vcf ];then
  if [ ! -s ${samfile%sam}sorted.vcf ]; then
    echo "Samtools: Creating pileup file from ${samfile%sam}sorted.bam"
    nice samtools mpileup -uf $reference_file ${samfile%sam}sorted.bam | bcftools view -vcg - > ${samfile%sam}sorted.vcf
  else
    echo "${samfile%sam}sorted.vcf exists. Remove or rename  to run pileup again. "
  fi
fi
}








function exit_on_error(){
  echo "Program $1 failed. Exiting Now";
  exit;

}
