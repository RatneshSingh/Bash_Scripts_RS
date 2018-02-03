### Remember to run S-7066 seperately as it was giving problem due to missing R2 read files.
### read file names and collect names of sample. Run trinity for each sample (leaf + stem combined)
function merge_files(){
    outfile=$1
    infiles=$2
    echo
    echo "***** MERGING FILES *****  "
    echo "output file: $outfile"
    echo "input files: $infiles"
    
    if  (echo $infiles|grep  -q ".gz$" ); then
        app="zcat"
    else
        app="cat"
    fi
    
    
    mkdir -p $wdir/${read_folder}_unified/$name
    
	if [[ -f ${outfile}.cat.working || -f ${outfile}.cat.OK ]]; then
            echo " "
            echo "Either Merged file exists or some other program is working on $outfile"
            echo "Skipping to next file"
            echo "**** Nothing to do for $outfile  ****"
            echo " "
            #rm $PWD/2_trimmed_reads_merged/$name/${tissue}${name}_All_fastq_merged.fastq.cat.working 
            #rm $PWD/2_trimmed_reads_merged/$name/${tissue}${name}_All_fastq_merged.fastq.cat.OK
     else
            echo ""
            if [[ -z $infiles ]];then echo; echo "No files detected in $tissue of $infiles samples ";continue; fi
            echo "Merging files $infiles  to $outfile "
            tempnum=$(date +%s).${RANDOM}
            touch ${outfile}.cat.working
            echo "running command:"
            echo "$app $infiles > $temp_folder/${tempnum}.fastq && sleep 1 &&  mv $temp_folder/${tempnum}.fastq $outfile && touch ${outfile}.cat.OK "
            $app $infiles > $temp_folder/${tempnum}.fastq && sleep 1 &&  mv $temp_folder/${tempnum}.fastq $outfile && touch ${outfile}.cat.OK 
            rm  ${outfile}.cat.working
            ####rm $PWD/2_trimmed_reads_merged/$name/${tissue}${name}_All_fastq_merged.fastq.cat.OK
            [[ -f $temp_folder/${tempnum}.fastq ]] && rm $temp_folder/${tempnum}.fastq	
     fi
    
    echo "***** END MERGING FILES *****"
    echo
}


########################################################################
##### 
estmeth=$1
: ${estmeth:="RSEM"}

cross_sample_norm=$2
: ${cross_sample_norm:="UpperQuartile"}   ### cross sample normalization methods e.g.  TMM, UpperQuartile, none

wdir="/mnt/nfs1/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations"
refname="LAP.US56.combined"  ### will be used to name folder and files.
transcript="$wdir/2_trimmed_reads_merged/${refname}/LAP.US56.9202.combined.fasta"

TRINITY_HOME="/usr/local/trinity"

read_folder="2_trimmed_reads"
#read_folder="1_raw_reads"

trinity="${TRINITY_HOME}/Trinity"
temp_folder="/home/ratnesh.singh/local_storage/tmp_folders"
local_wdir="/home/ratnesh.singh/local_storage"
matrix_fold="$wdir/7.1_Expression_Matrix"


### using trimmed paired reads
pt_command_list_file="$wdir/trinity_${read_folder}_on_${refname}_abundance_using_${estmeth}_per_tissue_commands_to_run.list"
ps_command_list_file="$wdir/trinity_${read_folder}_on_${refname}_abundance_using_${estmeth}_per_sample_commands_to_run.list"
pt_L_result_files_iso=()
pt_L_result_files_gene=()
pt_S_result_files_iso=()
pt_S_result_files_gene=()

ps_result_files_iso=()
ps_result_files_gene=()



rm $pt_command_list_file  $ps_command_list_file

#extra=""

### create index if it does not exists
if [[ ${estmeth,,} = "kallisto" &&  ! -f ${transcript}.kallisto_idx ]]     ;then 
     echo "**** CREATING KALLISTO INDEX ****"
     /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --prep_reference --est_method $estmeth $extra  --trinity_mode --output_dir $(dirname $transcript)
fi

if [[  ${estmeth,,} = "rsem"   &&    ! -f ${transcript}.RSEM.idx.fa  ]];then
    echo "**** CREATING RSEM INDEX ****"
    /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --prep_reference --est_method $estmeth $extra  --trinity_mode --aln_method bowtie2 --output_dir $(dirname $transcript)
fi  

if [[ $estmeth = "salmon"   && ! -f ${transcript}.salmon_quasi.idx ]];then
     echo "**** CREATING SALMON INDEX ****"
     /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --prep_reference --est_method $estmeth $extra  --trinity_mode --output_dir $(dirname $transcript)
fi

#### create bowtie2 index if alignment based methods are used for DEG
if [[ ${estmeth,,} = "rsem"   ||   ${estmeth,,} = "express"   ]] ;then
        extra="  --aln_method bowtie2   "
        if [[ ! -f ${transcript}.bowtie2.1.bt2  ]];then
	     echo "**** CREATING BOWTIE INDEX ****"
            /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --prep_reference --est_method $estmeth $extra  --trinity_mode  --aln_method bowtie2 --output_dir $(dirname $transcript)
        fi 
 fi
 



for name in $(ls $wdir/1_raw_reads/|grep "fastq"|cut -f1 -d "_"|cut -f2,3 -d "-"|sort|uniq|egrep   "(^7001|^7003|^7005|^7006|^7011|^7013|^7016|^7022|^7027|^7034|^7036|^7037|^7038|^7039|^7044|^7046|^7048|^7049|^7051|^7054|^7063|^7064|^7067|^7069|^7074|^7078|^7079|^7082|^7083|^7084|^7087|^7090|^7091|^7092|^7094|^7098|^7101|^7102|^7120|^7121|^7128|^9202|^LA-Purple|^US)"); do
    for tissue in "L-" "S-"; do
        ### used for untrimmed reads
        #out_fold="$wdir/7_estimate_abundance/$name/${tissue}${name}.trinity.abundance_originalUntrimmedReads_on_${refname}.$estmeth"
        
        ### used for trimmed paired reads
        out_fold="$wdir/7_estimate_abundance/$name/${tissue}${name}.trinity.abundance_${read_folder}_on_${refname}.$estmeth"
        out_name=${out_fold##*/}
	local_out_fold="$local_wdir/7_estimate_abundance/$name/${tissue}${name}.trinity.abundance_${read_folder}_on_${refname}.$estmeth"
    	mkdir -p $out_fold
	mkdir -p $local_out_fold
      	## Read all the file for sample $name
      	## for untrimmed reads
        #R1=$(ls $wdir/1_raw_reads/$tissue$name*R1*.gz);
        
        ## for trimmed paired reads
        R1=$(ls $wdir/${read_folder}/$name/$tissue$name*R1*P.qtrim.gz);
        
      	#R1_seqs=${R1[@]};
        #left=${R1_seqs// /,}
     	#R2_seqs=${R1_seqs//R1.fastq/R2.fastq}
    
        mkdir -p $wdir/${read_folder}_unified/$name    ### for trimmed reads
       ### merge reads for R1 and R2 seperately to combine multiple sequencing results for each tissue type.
      for reads in R1; do
            treads_R=$(ls $wdir/${read_folder}/$name/$tissue$name*${reads}*.P.qtrim.gz);
            treads_R1=${treads_R//$'\n'/  }
            treads_R2=${treads_R1//R1/R2}
            merge_files  "$wdir/${read_folder}_unified/$name/${tissue}${name}_R1_all.fastq" "$treads_R1"
            merge_files  "$wdir/${read_folder}_unified/$name/${tissue}${name}_R2_all.fastq" "$treads_R2"
           #merge_files  "$wdir/${read_folder}_unified/$name/${tissue}${name}_${reads}_all.fastq" "${reads_R//$'\n'/  }"
           #merge_files  "$wdir/${read_folder}_unified/$name/${tissue}${name}_${reads}_all.fastq" "${reads_R//$'\n'/  }"
       done
        ### create commandlist to be run in parrallel
         if [[ -f  $wdir/${read_folder}_unified/$name/${tissue}${name}_R1_all.fastq && -f  $wdir/${read_folder}_unified/$name/${tissue}${name}_R2_all.fastq ]];then
	     if [[ ! -f ${out_fold}/${out_fold##*/}.isoforms.results.ok ]]; then               
		echo   "nice /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --est_method $estmeth   --trinity_mode --seqType fq --thread_count 30 --output_dir $local_out_fold  --output_prefix ${out_name}  --left  $wdir/${read_folder}_unified/$name/${tissue}${name}_R1_all.fastq  --right $wdir/${read_folder}_unified/$name/${tissue}${name}_R2_all.fastq  $extra > ${out_fold}.log  2>&1 && cp -fr $local_out_fold ${out_fold%/*}/ " | tee --append   $pt_command_list_file
                
	     else
		echo "***** The following output result file exists. skipping this sample:"
		echo ${out_fold}/${out_fold##*/}.isoforms.results.ok	
	     fi
     else
              echo " One or both of the Left and right file does not exists"
               ls -h $wdir/${read_folder}_unified/$name/${tissue}${name}_R1_all.fastq $wdir/${read_folder}_unified/$name/${tissue}${name}_R2_all.fastq
         fi
        
    
    done
    
    
    ## merge for combined run per sample instead of  per tissue
    ## merge_files   cat|zcat     "outputfile_wid_full_path"   "reads"
    for reads in R1 ; do
         preads_R=$(ls $wdir/${read_folder}/$name/*${name}*${reads}*.P.qtrim.gz)
         preads_R1=${preads_R//$'\n'/  }
         preads_R2=${preads_R1//R1/R2}
         merge_files  "$wdir/${read_folder}_unified/$name/${name}_R1_all.fastq" "$preads_R1"
         merge_files  "$wdir/${read_folder}_unified/$name/${name}_R2_all.fastq" "$preads_R2"
    done
        
    ### create commandlist to be run  using parrallel
    out_fold="$wdir/7_estimate_abundance/$name/${name}.trinity.abundance_${read_folder}_on_${refname}.$estmeth"
    out_name=${out_fold##*/}
    local_out_fold="$local_wdir/7_estimate_abundance/$name/${name}.trinity.abundance_${read_folder}_on_${refname}.$estmeth"
    
    if [[ -f  $wdir/${read_folder}_unified/$name/${name}_R1_all.fastq && -f  $wdir/${read_folder}_unified/$name/${name}_R2_all.fastq ]];then
	if [[ ! -f ${out_fold}/${out_fold##*/}.isoforms.results.ok ]]; then 
        	echo   "nice /usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $transcript --est_method $estmeth --trinity_mode --seqType fq --thread_count 30 --output_dir $local_out_fold  --output_prefix ${out_name}  --left  $wdir/${read_folder}_unified/$name/${name}_R1_all.fastq  --right $wdir/${read_folder}_unified/$name/${name}_R2_all.fastq  $extra > ${out_fold}.log   2>&1 && cp -fr $local_out_fold ${out_fold%/*}/" | tee --append   $ps_command_list_file
        
	else
		echo "***** The following output result file exists. skipping this sample:"
		echo ${out_fold}/${out_fold##*/}.isoforms.results.ok	
	fi
    else
           echo " One or both of the Left and right file does not exists"
           ls -h $wdir/${read_folder}_unified/$name/${name}_R1_all.fastq $wdir/${read_folder}_unified/$name/${name}_R2_all.fastq
    fi
    
     
done #> trinity_merged_abundance_per_tissue_commands_to_run.list 




## run te command in parrallel
#
#rm All_commands_to_run.list ; cat $ps_command_list_file $pt_command_list_file > All_commands_to_run.list
#parallel  -j 10 -a $ps_command_list_file 
parallel  -j 10 -a $pt_command_list_file
#parallel  -j 10 < All_commands_to_run.list
#echo "**** EXITING AFTER COMMAND CREATION **** ";echo "**** UNCOMMENT LINE $LINENO IN " $(basename $0) " TO RUN MATRIX CREATION STEP ****";exit   ## uncomment to stop it from running further

#exit;

### create a expression matrix from all the sugarcane samples
printf "\n\n**** CREATING EXPRESSION MATRIX ****\n"
if [[ ${estmeth,,} = "kallisto" ]];then
generesfile="abundance.tsv.genes"
isoresfile="abundance.tsv"
elif [[ ${estmeth,,} = "rsem" ]];then
generesfile="genes.results"
isoresfile="isoforms.results"
elif [[ ${estmeth,,} = "salmon" ]];then
generesfile="quant.sf.genes"
isoresfile="quant.sf"
elif [[ ${estmeth,,} = "express" ]];then
generesfile="results.xprs.genes"
isoresfile="results.xprs"
fi


out_fold="$wdir/7_estimate_abundance/$name/${tissue}${name}.trinity.abundance_${read_folder}_on_${refname}.$estmeth"
out_name=${out_fold##*/}
reserve_generesult_for_seqlength="None"
reserve_isoresult_for_seqlength="None"
#matrix_fold="$wdir/7.1_Expression_Matrix"  ## moved it up
mkdir -p $matrix_fold
for restype in  $generesfile  $isoresfile ; do
    L_list=""
    S_list=""
    All_list=""
    mtype=""

    L_list=$(ls $wdir/7_estimate_abundance/*/L-*trinity.abundance_${read_folder}_on_${refname}.$estmeth/*${restype})
    S_list=$(ls $wdir/7_estimate_abundance/*/S-*trinity.abundance_${read_folder}_on_${refname}.$estmeth/*${restype})
    All_list=$(ls $wdir/7_estimate_abundance/*/*trinity.abundance_${read_folder}_on_${refname}.$estmeth/*${restype}|egrep -v "(L-|S-)" )
 	if ($(echo "$restype" |grep -q "gene"));then mtype="genes";else mtype="isoform"; fi;echo $mtype
  if ($(echo "$restype" |grep -q "gene"));then reserve_generesult_for_seqlength=(${L_list//$'\n'/  }); fi
 if ($(echo "$restype" |grep -q "isoform"));then reserve_isoresult_for_seqlength=(${L_list//$'\n'/  }); fi
	echo "processing results from $estmeth for $mtype "
#### matrix creation for leaf samples.
    if [[ ! -f   $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Leaf.${mtype}.${cross_sample_norm}.EXPR.matrix ]]; then    /usr/local/trinity/util/abundance_estimates_to_matrix.pl  --est_method $estmeth  --cross_sample_norm  $cross_sample_norm  --name_sample_by_basedir  --out_prefix  $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Leaf.${mtype}  ${L_list//$'\n'/  }; 


else
	echo " **** MATRIX FILE: "$matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Leaf.${mtype}.${cross_sample_norm}.EXPR.matrix " exists."
	echo "*** SKIPPING MATRIX CREATION STEP ****"; echo
fi


#### matrix creation for stem samples.
   if [[ ! -f   $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Stem.${mtype}.${cross_sample_norm}.EXPR.matrix ]] ;then  /usr/local/trinity/util/abundance_estimates_to_matrix.pl  --est_method $estmeth  --cross_sample_norm  $cross_sample_norm  --name_sample_by_basedir  --out_prefix  $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Stem.${mtype}   ${S_list//$'\n'/  }
else
	echo " **** MATRIX FILE: " $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_Stem.${mtype}.${cross_sample_norm}.EXPR.matrix " exists."
	echo "*** SKIPPING MATRIX CREATION STEP ****"; echo
fi



#### matrix creation for LeafStem samples.
   if [[ ! -f   $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LeafStem.${mtype}.${cross_sample_norm}.EXPR.matrix ]] ;then  /usr/local/trinity/util/abundance_estimates_to_matrix.pl  --est_method $estmeth  --cross_sample_norm  $cross_sample_norm  --name_sample_by_basedir  --out_prefix  $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LeafStem.${mtype}   ${S_list//$'\n'/  } ${L_list//$'\n'/  }
else
	echo " **** MATRIX FILE: " $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LeafStem.${mtype}.${cross_sample_norm}.EXPR.matrix " exists."
	echo "*** SKIPPING MATRIX CREATION STEP ****"; echo
fi






#### matrix creation for leaf+stem merged samples.
   if [[ ! -f   $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LSCombined.${mtype}.${cross_sample_norm}.EXPR.matrix ]]  ;then   /usr/local/trinity/util/abundance_estimates_to_matrix.pl  --est_method $estmeth  --cross_sample_norm  $cross_sample_norm  --name_sample_by_basedir  --out_prefix  $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LSCombined.${mtype}   ${All_list//$'\n'/  } 
else
	echo " **** MATRIX FILE: "$matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_LSCombined.${mtype}.${cross_sample_norm}.EXPR.matrix " exists."
	echo "*** SKIPPING MATRIX CREATION STEP ****"; echo
fi






### Shorten the name of header columns. 
	for mfile in $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_*.${mtype}*; do 
		if $(grep -iq ".trinity.abundance_${read_folder}_on_${refname}.${estmeth}" $mfile) ; then sed -i "s/.trinity.abundance_${read_folder}_on_${refname}.${estmeth}//g" $mfile ; fi
		if $(grep -iq ".trinity.abundance_originalUntrimmedReads_on_${refname}.${estmeth}" $mfile)  ; then sed -i "s/.trinity.abundance_originalUntrimmedReads_on_${refname}.${estmeth}//g" $mfile; fi
	done


	#echo  ${L_list//$'\n'/ }
	#echo  ${S_list//$'\n'/ }
	#echo  ${All_list//$'\n'/ } 


done

echo "**** CREATING ExN50 TABLE ****"
### estimating ExN50 for the assembled genes.

for restype in  $generesfile  $isoresfile ; do
	if ($(echo "$restype" |grep -q "gene"));then mtype="genes";else mtype="isoform"; fi;echo $mtype

	for mfile in $matrix_fold/${estmeth}_express.matrix_${read_folder}_On${refname}_NormBy.${cross_sample_norm}_*.${mtype}.${cross_sample_norm}.EXPR.matrix; do 
		if [[ ! -f $mfile.ExN50.stats ]] ; then 
			$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl  $mfile $transcript | tee $mfile.ExN50.stats 
		else 
			echo "**** SKIPPING ExN50 RUN. OUTPUT EXISTS :$mfile.ExN50.stats .  ****" 
		fi
	done
done


##### create fasta file of longest isoform
[[ -e ${transcript%.*}.longest_isoform.fasta ]] || (/usr/local/trinity/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $transcript > ${transcript%.*}.longest_isoform.fasta && sed -r "s/_i\d*.*//g" ${transcript%.*}.longest_isoform.fasta > ${transcript%.*}.longest_isoform_as_genes.fasta)



#### create seq_lengths table for GO annotation step in (_DiffExp_analysis.sh
if ( -e ${transcript}.gene_iso.seq_lengths ); then echo "${transcript}.gene_iso.seq_lengths exists"; 
else 
	cat ${reserve_isoresult_for_seqlength[0]} | cut -f 1,3  > ${transcript}.gene_iso.seq_lengths
	cat ${reserve_generesult_for_seqlength[0]}|tail -n +2 >> ${transcript}.gene_iso.seq_lengths
	echo "Gene lengths saved in ${transcript}.gene_iso.seq_lengths"
fi
#echo ${reserve_isoresult_for_seqlength[0]} ${reserve_generesult_for_seqlength[0]}|cut -f 1,3  > ${transcript}.gene_iso.seq_lengths
#cat ${reserve_isoresult_for_seqlength[0]} ${reserve_generesult_for_seqlength[0]}|cut -f 1,3  > ${transcript}.gene_iso.seq_lengths

