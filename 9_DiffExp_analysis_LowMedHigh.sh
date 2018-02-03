echo;echo;
echo "******************************************************************************"
echo " **** START OF " $0 " SCRIPT *****"
echo "******************************************************************************"

##### values needed from commandline
matrix=$1    # count matrix
group_list=$2
demeth=$3
cross_sample_norm=$4
#ref_sample=$5
#annot=$4


group_list=$(realpath $group_list)
if [[ -z $matrix  ||  -z $group_list ]];then echo "usage: script  raw_read_count_matrix group_list  est_meth[DESeq2|edgeR|voom|ROTS]";  exit ; fi

if [[ ! -f $matrix  ||  ! -f $group_list ]];then echo "Unable to open matrix or group file:";echo "usage: script  count_matrix group_list  est_meth[DESeq2|edgeR|voom|ROTS]";  exit ; fi


### Some variables required for the analysis.
refname="LAP.US56.combined"  ### will be used to name folder and files.
#cross_sample_norm="UpperQuartile"
feat_map="8_annotate_assemblies/${refname}_trimmed_merged_assembly/annot_feature_map.txt"
wdir="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations"
diff_exp_folder="9_DiffExp_analysis_2trimOnLAP.US.9202CombinedAssembly_LowHighsmallerGrpFinal"
 ################################################################################################################
#####  MOST PROBABLY NO CHANGE NEEDED BELOW THIS
##########
TRINITY_HOME="/usr/local/trinity"
: ${demeth:="DESeq2"}
: ${cross_sample_norm:="TMM"}
qval=0.05
log2foldc=1

matrix_file=$(basename $matrix)
matrix_path=$(realpath $(dirname $matrix))
group=$(basename $group_list)
diffold=${demeth}.wAnnot.${matrix_file%.*}.grp.${group%.*}
script_folder="/home/ratnesh.singh/Scripts/Perl"

### to delete existing run.
#rm -r $wdir/${diff_exp_folder}/$diffold
###
[[ -d $wdir/${diff_exp_folder} ]] || mkdir $wdir/${diff_exp_folder}

outprefix="$wdir/${diff_exp_folder}/$diffold/${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix}.diffExpr.qval${qval}_Log2FC${log2foldc}"
transcript="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/8_annotate_assemblies/${refname}_trimmed_merged_assembly/LAP.US56.9202.combined.fasta"
goannot="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/8_annotate_assemblies/${refname}_trimmed_merged_assembly/wAnnot.go_annotations.txt"
genegoannot="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/8_annotate_assemblies/${refname}_trimmed_merged_assembly/go_annotations_genes.txt"
gene_length=$(dirname $transcript)/wAnnot.$(basename $transcript).gene_iso.seq_lengths
#gene_length="$wdir/wAnnot.LAP.US56.9202.combined.fasta.regions.All.firsthead_streamsort.ReadCount_DPgt500.OneLinePerAllele.counts.matrix.All_seq_lengths_as_1"

longest_isoform_as_genes="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/2_trimmed_reads_merged/${refname}/LAP.US56.9202.combined.longest_isoform_as_genes.fasta"
## if gene length does not exist, create one
[[ -e $(dirname $transcript)/wAnnot.$(basename $transcript).gene_iso.seq_lengths ]] || $TRINITY_HOME/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl ${transcript}.gene_iso.seq_lengths $(dirname $transcript)/annot_feature_map.txt > $(dirname $transcript)/wAnnot.$(basename $transcript).gene_iso.seq_lengths

mkdir -p $wdir/${diff_exp_folder}/$diffold
cp $group_list $wdir/${diff_exp_folder}/$diffold/



tree_clust_percent=60
DE_for=$( echo $(sed -r 's/\s+/\t/g' $group_list|cut -f 1|sort|uniq)|sed -r 's/\s+/_vs_/g')
samp_names=($( echo $(sed -r 's/\s+/\t/g' $group_list|cut -f 1|uniq)|sed -r 's/\s+/ /g'))
#: {ref_sample:=${samp_names[1]}}
#ref_sample=${samp_names[1]}

ref_sample=$(echo ${group/.txt/}|cut -d "_" -f 4)


contrast_file="$wdir/${diff_exp_folder}/${diffold}.contrast.txt"
contrast=($(echo ${group/.txt/}|cut -d "_" -f 4,6|sed -r 's/_/\t/g'))
printf "${contrast[0]}\t${contrast[1]}" > $contrast_file


#exit
##
## Annotate all matrix before running DE analysis so annotation could be exported with the analysis result.

for mat_file in $matrix_path/${matrix_file/.counts.matrix/}*; do
[[ -S $matrix_path/wAnnot.$(basename $mat_file) ]] || $TRINITY_HOME/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl $mat_file  $feat_map > $matrix_path/wAnnot.$(basename $mat_file)
done

## run DESeq analysis on matrix containing RAW read counts.
echo "**** RUNNING DE ANALYSIS AT LINE: $LINENO ****"
if [[ ! $(ls $wdir/${diff_exp_folder}/$diffold/wAnnot.${matrix_file}.*.${demeth}.DE_results ) ]] ; then
[[ -s $wdir/${diff_exp_folder}/$diffold/wAnnot.${matrix_file}.*.${demeth}.DE_results ]] || $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ${matrix_path}/wAnnot.${matrix_file} --method $demeth --samples_file $group_list --contrast $contrast_file  --output $wdir/${diff_exp_folder}/$diffold




else
 echo " output directory $wdir/${diff_exp_folder}/$diffold exists ";
 echo " Skipping run_DE_analysis step."
 echo "following is the content of output folder. Delete the folder to run again."
 ls -l $wdir/${diff_exp_folder}/$diffold

echo;echo

fi








#exit;
## extract and analyze diff exp genes
curdir=$(pwd)
echo "**** ANALYZING DIFFERENTIALLY EXPRESSED GENES AT $LINENO ****"
cd $wdir/${diff_exp_folder}/$diffold || (echo "Unable to find $wdir/${diff_exp_folder}/$diffold ";exit)

#### if TMM normalized matrix file does not exists... Create one.
if [[ ! -e ${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix} ]];then
	perl /usr/local/trinity/util/support_scripts/run_TMM_scale_matrix.pl --matrix ${matrix_path}/wAnnot.${matrix_file} > ${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix}
#else
#  echo "Skipping: ${cross_sample_norm} Normalized expression matrix exists:${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix"

fi

if [[ ! -e ${outprefix}.matrix.log2.dat ]] ; then
	if [[ -e ${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix} ]];then
		$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix} --output $outprefix --samples $group_list --examine_GO_enrichment --GO_annots $goannot  --gene_lengths $gene_length -P $qval -C $log2foldc
	else
		echo "Unable to run analyze_diff_expr.pl at $LINENO";echo "Normalized TMM Matrix file does not exists: " ${matrix_path}/wAnnot.${matrix_file/counts.matrix/${cross_sample_norm}.EXPR.matrix}

	fi
else
	echo " **** DIFFERENTIAL EXPRESSION RESULT FILE FROM PREVIOUS RUN: "${outprefix}.matrix.log2.dat " exists."
	echo "*** SKIPPING DIFFERENTIAL EXPRESSION ANALYSIS STEP ****"; echo
fi

###########################################################################################################
#### create fold change table from DE genes.
echo "**** CREATING FOLD CHANGE TABLE AND FASTA FILE AT LINE: $LINENO ****"
for subset in *UP.subset; do
    [[ -e $subset ]] || continue
    if [[ ! -s ${subset/UP.subset/UP.FC.table} ]] ; then
	cut -f 1,2,3,4,5 $subset > ${subset/UP.subset/UP.FC.table}
    else
	echo Fold change table exists. Skipping : ${subset/UP.subset/UP.subset.FC.table}
    fi


    ### create fasta file of up regulated sequences
    if [[ ! -s ${subset/UP.subset/UP.FC.fasta} ]] ; then
	(cut -f 1 $subset|sed 's/\^.*//' > ${subset/UP.subset/UP.FC.list}) && (perl $script_folder/search_seq-V3.4.pl -l ${subset/UP.subset/UP.FC.list} -o ${subset/UP.subset/UP.FC.fasta} -s $longest_isoform_as_genes -m exact -d " " -c 1)
    else
	echo Fasta sequence file exists. Skipping : ${subset/UP.subset/UP.FC.fasta}
  fi

######################################################
  #### create a mapping table for DE genes with RefSeq or other databases. Some enrichment analysis programs need these names.
  #sed -r "s/\^\S+//" ${subset/UP.subset/UP.FC.table} > ${subset/UP.subset/UP.FC.forMapByBlast.table}
   perl $script_folder/Map_sequence_by_blast.pl -qs ${subset/UP.subset/UP.FC.fasta} -qt ${subset/UP.subset/UP.FC.table} -ts /mnt/nfs1/Databases/BlastDatabases/NCBI/refseq_protein -ts /mnt/nfs1/Databases/BlastDatabases/NCBI/swissprot -ts /mnt/nfs1/Databases/BlastDatabases/NCBI/pdbaa -diamond -bt blastx -delim '^' -col 1

###nn
done


#### create table for mapman.
echo "**** CREATING TABLE FOR MAPMAN ****"
for de_result in wAnnot.${matrix_file}.*.${demeth}.DE_results; do
    if [[ -e $de_result && ! -s ${de_result/DE_results/}forMapMan.qval0.05.txt ]] ; then
	    if [[ $(head -n 1 $de_result|cut -f1) -eq 'id' ]]; then
        $(echo $demeth|grep -q -i edgeR) || cut -f 1,5,8,9 ${de_result} |sed -r "s/\^\S+//g" > ${de_result/DE_results/}forMapMan.txt
        $(echo $demeth|grep -q -i edgeR) && sed -r "s/\s+/\t/g" ${de_result}|sed -r "s/^logFC/id\tlogFC/g"|cut -f 1,2,4,5|sed -r  "s/\^\S+//g"  > ${de_result/DE_results/}forMapMan.txt
		else
			### the newer version of trinity has different number of columns. modified lines below to needed info.
			$(echo $demeth|grep -q -i edgeR) || head -n 1 ${de_result}|cut -f 6,9,10 |sed -r "s/^\s*/id\t/" > ${de_result/DE_results/}forMapMan.txt
			$(echo $demeth|grep -q -i edgeR) || tail -n +2 ${de_result}|cut -f 1,7,10,11 |sed -r "s/\^\S+//g" >> ${de_result/DE_results/}forMapMan.txt
		fi
	$(echo $demeth|grep -q -i edgeR) && sed -r "s/\s+/\t/g" ${de_result}|sed -r "s/^logFC/id\tlogFC/g"|cut -f 1,2,4,5|sed -r  "s/\^\S+//g"  > ${de_result/DE_results/}forMapMan.txt
  ### create seperate tables for qvalues
	head -n 1 ${de_result/DE_results/}forMapMan.txt > ${de_result/DE_results/}forMapMan.qval0.05.txt;  awk '{if($4 <= 0.05) print }' ${de_result/DE_results/}forMapMan.txt |sed -r "s/\^\S+//g" >> ${de_result/DE_results/}forMapMan.qval0.05.txt
	head -n 1 ${de_result/DE_results/}forMapMan.txt > ${de_result/DE_results/}forMapMan.pval0.05.txt;  awk '{if($3 <= 0.05) print }' ${de_result/DE_results/}forMapMan.txt |sed -r "s/\^\S+//g" >> ${de_result/DE_results/}forMapMan.pval0.05.txt

    else
	echo MapMan table exists. Skipping : ${de_result/DE_results/}forMapMan.txt
    fi
done

#### create table of gene names go terms for input for wego
for de_result in wAnnot.${matrix_file}.*.${demeth}.DE_results; do
    if [[ -e $de_result && ! -s ${de_result/DE_results/}forWeGo.qval0.05.txt ]] ; then

	cut -f1 ${de_result/DE_results/}forMapMan.qval0.05.txt > ${de_result/DE_results/}forWeGo.qval0.05.list
	perl -e '
		open $data,$ARGV[0];my%godata;
		while(<$data>){
		    my($gname,$rest)=split /\s+/;
		    $godata{$gname}=$_;

		}

		close $data;
		open $list,$ARGV[1];
		open $out,">",$ARGV[2];
		while (my$lname = <$list>){
		    $lname=~s/\s+//g;
		    print $out $godata{$lname};

		}
		close $out;close $list;


	' $genegoannot ${de_result/DE_results/}forWeGo.qval0.05.list ${de_result/DE_results/}forWeGo.qval0.05.txt
 	#
	#grep -f ${de_result/DE_results/}forWeGo.qval0.05.list $genegoannot > ${de_result/DE_results/}forWeGo.qval0.05.txt
	sed -i -r "s/,/\t/g" ${de_result/DE_results/}forWeGo.qval0.05.txt
    fi

   if [[ -e ${de_result/DE_results/}forWeGo.qval0.05.txt && ! -s ${de_result/DE_results/}forWeGo.qval0.05.table ]]; then
	perl -e 'while(<>){s/^\s+//;my($gene,@goid)=split /\s+/; foreach my$gid(@goid){print "$gene\t$gid\n"}}' < ${de_result/DE_results/}forWeGo.qval0.05.txt > ${de_result/DE_results/}forWeGo.qval0.05.table
   fi


done






#### creating Summary report:
## create summary table of up and down regulated genes from *.DE.results
 #sum_file="Summary5.txt.table"

for de_result in wAnnot.${matrix_file}.*.${demeth}.DE_results; do
  if [[ -e $de_result && ! -s ${de_result/DE_results/}DE_summary.txt ]] ; then
  	res=${de_result}
  	if [[ -e ${de_result} ]]; then
       perl ~/Scripts/Perl/process_DE_results_summary.pl  ${de_result} ${de_result/DE_results/}DE_summary.txt

		     ####
		       sed -i -r "s/File/Tissue\tComparison/g" ${de_result/DE_results/}DE_summary.txt
		       sed -i -r "s/\S+\///g" ${de_result/DE_results/}DE_summary.txt
		       sed -i -r "s/wAnnot.RSEM_express.matrix_2_trimmed_reads_On9202_NormBy.TMM_//g" ${de_result/DE_results/}DE_summary.txt
		       sed -i  -r "s/.genes.counts.matrix./\t/g" ${de_result/DE_results/}DE_summary.txt
		       sed -i  -r "s/.DE_results//g" ${de_result/DE_results/}DE_summary.txt
    fi
  fi
done




#### create list of enriched go terms filtered by FDR.
for subset in *subset.GOseq.enriched; do
  [[ -z $subset ]] && continue
  awk 'FS="\t" { if(NR == 1 || $8 <= 0.05) print } '  $subset |cut -f1 > ${subset}.FDR_0.05.list;
  awk 'FS="\t" { if(NR == 1 || $8 <= 0.05) print } '  $subset > ${subset}.FDR_0.05.table;
done

#### seperate enriched GO terms based on FDR.
for subset in *.subset.GOseq.depleted; do
	[[ -z $subset ]] && continue
  	awk 'FS="\t" { if(NR == 1 || $9 <= 0.05) print } '  $subset |cut -f1 > ${subset}.FDR_0.05.list;
	awk 'FS="\t" { if(NR == 1 || $9 <= 0.05) print } '  $subset > ${subset}.FDR_0.05.table;
done







###########################################################################################################
### cut trees in groups
echo "**** Not RUNNING CUT GROUPS FROM TREE AT $LINENO IN " $(basename $0) " ****";
#[[ -d ${outprefix}.matrix.RData.clusters_fixed_P_$tree_clust_percent ]] || $TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R  ${outprefix}.matrix.RData --Ptree $tree_clust_percent















#echo "**** SKIPPING KO ANNOTATION STEP ****"; echo "**** COMMENT $LINENO IN " $(basename $0) " TO RUN KAAS ANALYSIS ON WEB SERVER **** ";exit;
### get KO annotation for enriched genes

##create a list of names for upregulated genes

echo
echo "**** EXTRACTING FASTA SEQUENCES FOR KO ANNOTATION **** "
for subset in $wdir/${diff_exp_folder}/$diffold/wAnnot.${matrix_file}.*.${demeth}.DE_results.P${qval}_C${log2foldc}.*.subset; do
	cut -f1 $subset|cut -f1 -d "^" > ${subset}.list
	if $(echo ${matrix_file}|grep -q genes );then
		if [[ ! -f ${subset}.fasta ]]; then
			perl ~/Scripts/Perl/search_seq-V3.4.pl -s $longest_isoform_as_genes -l ${subset}.list -o ${subset}.fasta -d "\s+" -c 1  ## original
			#perl ~/Scripts/Perl/search_seq-V3.4.pl -s $longest_isoform_as_genes -l ${subset}.list -o ${subset}.fasta -d ":" -c 1     ## modified for ASE
		else
			echo "*** SUBSET FASTA FILE EXISTS:"${subset}.fasta "SKIPPING SEQUENCE EXTRACTION STEP ****"
		fi
	else
		if [[ -f ${subset}.fasta ]] ; then
			perl ~/Scripts/Perl/search_seq-V3.4.pl -s ${transcript} -l ${subset}.list -o ${subset}.fasta -d " " -c 1  ## original
			#perl ~/Scripts/Perl/search_seq-V3.4.pl -s ${transcript}  -l ${subset}.list -o ${subset}.fasta -d ":" -c 1     ## modified for ASE
		else
			echo "*** SUBSET FASTA FILE EXISTS:"${subset}.fasta "SKIPPING SEQUENCE EXTRACTION STEP ****"
		fi
	fi
done

echo
echo "**** RUNNING KAAS ANNOTATION STEP **** "
#sh ~/Scripts/Bash_scripts/get_png_kaas.sh  $(ls *high_biomass-UP.subset.fasta) $(ls *low_biomass-UP.subset.fasta)
echo;echo

#########################################################################################################################
##### finding Allele specific expression for DE genes.
cd $curdir
ploidy=12
min_coverage=3
picard="java -jar /usr/local/picard/picard.jar "
#### if provided group list do for group only.
if [[ ! -z $group_list ]]; then
    grp=$group_list
    echo $grp
    grp_list=$(cut -f2 $grp|sort|uniq|grep -v -e "  ")
    #echo $grp_list
    grp_list=${grp_list//$'\n'/ }
    #grp_list=${grp_list//$'\n'/}
    #echo $grp_list
else
    echo "Usage $0 group_list";
    echo "Group list_not_provided. Using all the samples";
   grp_list=$(ls $wdir/7_estimate_abundance/)
fi




ase_fold="$wdir/${diff_exp_folder}/$diffold/ase_analysis"
#de_fold=$ase_fold/DE_analysis
###rm -r $ase_fold
mkdir -p $ase_fold
#mkdir -p $de_fold
##map reads on DE genes.
for subset in $(ls $wdir/${diff_exp_folder}/$diffold/wAnnot.${matrix_file}.*.${demeth}.DE_results.P${qval}_C${log2foldc}.*.subset.fasta); do
	#mat_file_base=$(basename $mat_file)
	#mat_file_base=${mat_file_base%%TPM*/}


	de_fold=${subset%%.subset.fasta}
	de_fold=${de_fold##*DE_results.P${qval}_C${log2foldc}.}
	de_fold=$ase_fold/${de_fold/-/_}_DE
	#exit
	short_outfile=$(echo $(basename $subset)|sed "s/wAnnot.${matrix_file%.*}.matrix.//g")
	### create index if it does not exists
		if [[  ${estmeth,,} = "rsem"   &&    ! -f ${subset}.RSEM.idx.fa  ]];then
   			 echo "**** CREATING RSEM INDEX ****"
    			/usr/local/trinity/util/align_and_estimate_abundance.pl --transcripts $subset --prep_reference --est_method $estmeth $extra  --trinity_mode --aln_method bowtie2
		fi

		#### create bowtie2 index if alignment based methods are used for DEG
	        if [[ ! -f $ase_fold/$(basename $subset).bowtie2.1.bt2  ]];then
			     echo "**** CREATING BOWTIE INDEX ****"
		            bowtie2-build $subset $ase_fold/$(basename $subset).bowtie2
		else
			echo "Bowtie2 index for $subset exists"
		fi

	echo "Running bowtie2 on $grp_list"
	for sample in $grp_list S-LA-Purple L-LA-Purple S-9202 L-9202 ;do
		### align reads to DE seq using bowtie2
		if [[ ! -e $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam ]]; then
			nice bowtie2 -x $ase_fold/$(basename $subset).bowtie2 -1 $(echo $wdir/2_trimmed_reads/${sample#*-}/${sample}*R1.fastq.gz.P.qtrim.gz|tr " " ,) -2 $(echo $wdir/2_trimmed_reads/${sample#*-}/${sample}*R2.fastq.gz.P.qtrim.gz|tr " " ,) -U $(echo $wdir/2_trimmed_reads/${sample#*-}/${sample}*R1.fastq.gz.U.qtrim.gz|tr " " ,)  -U $(echo $wdir/2_trimmed_reads/${sample#*-}/${sample}*R2.fastq.gz.U.qtrim.gz|tr " " ,)  --rg-id $sample  --rg LB:lib_${sample} --rg PL:illumina --rg PU:unit${sample} --rg SM:${sample}    --threads 48 --mm --qc-filter -S $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sam --no-unal
		fi

		### echo sort sam and save as bam.
		if [[ ! -e $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam ]]; then
			nice /usr/local/bin/samtools sort --threads 20 -m 10G -o $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam  $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sam && rm $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sam
		else
			echo "Bam alignment exists"
		fi


		#echo "Adding Readgroup to bam alignment"
		#if [[ ! -e $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam ]]; then
		#	$picard AddOrReplaceReadGroups I=$ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.bam  O=$ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam RGID=${sample} RGLB=lib_${sample} RGPL=illumina RGPU=unit${sample} RGSM=${sample}
			#ln -s $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.bam $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam
		#fi



		echo "indexing bam"
		if [[ ! -e  $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam.bai ]]; then
			echo "Indexing bam file"
			nice /usr/local/bin/samtools index  $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam
		else
			echo "Bam index for $sample exists"
		fi



	done



	echo "finding SNP using freebayes"
		### create bam list
		bam_list=" "
		#for sample in $grp_list;do
		#	bam_list="$bam_list $ase_fold/$(basename $subset)_vs_${sample}.bowtie2.sorted.rg.bam"
		#done

        echo "Using all the bam files in folder $ase_fold. Modify around here to use grp list"
 		bam_list=" "
		for sample in `ls $ase_fold |grep bowtie2.sorted.rg.bam$|grep ${subset##*/}`;do
			bam_list="$bam_list $ase_fold/$sample"
            echo "added $sample in list"
		done

#exit
		if [[ ! -s  $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf ]]; then
			nice freebayes --fasta-reference $subset --ploidy $ploidy --use-best-n-alleles 4 --pooled-continuous --min-coverage $min_coverage -F 0.1 --vcf $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf $bam_list

		else
			echo "VCF output $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf exists"

		fi


	echo "Process VCF to get read count per SNP"

	if [[ ! -s $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf.ReadCount.table ]]; then
		bcftools query -H -f '%CHROM\:%POS\:%REF\t[\t%RO]\n%CHROM\:%POS\:%ALT\t[\t%AO]\n' -o $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf.ReadCount.table $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf
	else
		echo "read count table exists:"

	fi


	### create Non-normalized ReadCount table one allele per line.
	if [[ ! -s $ase_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix ]]; then
		##create one line per allele | remove number from header | remove redundant header line > save to New file.

		perl -e 'while(<>){ if ($_ !~ m/,/){s/\./0/g;print "$_"; next} else{my$line=$_;my@elm=split /\s+/;my@id=split /:/,$elm[0];my@gts=split /,/,$id[-1];my$ntitle=join ":",@id[0 .. $#id-1];	for(my$i=0;$i<=$#gts;$i++){print "$ntitle:$gts[$i]";for(my$j=1;$j<=$#elm;$j++){my@cur_elm=split /,/,$elm[$j];my$rc = $cur_elm[$i] ? $cur_elm[$i] : 0;$rc= 0 if ($rc =~m/^\.$/ );print "\t$rc";} print "\n";}}}' < $ase_fold/$(basename $subset).bowtie2.sorted.rg.bam.freebayes.vcf.ReadCount.table |sed -r 's/\[[0-9]*\]//g'| sed -r '/^\s*CHROM:POS:ALT/d'|sed 's/:RO//g'|sed -r 's/^#\s*CHROM:POS:REF\s*/\t/'|sed -r 's/\t+/\t/g' > $ase_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix


	else
		echo "$ase_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix exists"
	fi



	### Create Normalized Read count TMM
	if [[ ! -s $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix ]]; then
		perl /usr/local/trinity/util/support_scripts/run_TMM_scale_matrix_RSmod.pl --matrix ${matrix_path}/wAnnot.${matrix_file} --select $ase_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix > ${de_fold}/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix

	else
		echo "Normalized Read counts exists: $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix"
	fi

 	### run DE expression analysis for normalized AS expression matrix.

	if [[ ! $(ls $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix.*.DE_results) ]] ; then
		$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix $ase_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix --method $demeth --samples_file $group_list --contrast $contrast_file  --output $de_fold
	fi


	### process DE result
	for de_results in $(ls $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.matrix.*.DE_results);do

		[[ -z ${de_results+x} ]] || cut -f 1,5,8,9 $de_results |sed -r "s/\^\S+//g" > ${de_results/DE_results/}forMapMan.txt
		[[ -z ${de_results+x} ]] || head -n 1 ${de_results/DE_results/}forMapMan.txt > ${de_results/DE_results/}forMapMan.qval0.05.txt;  awk '{if($4 <= 0.05) print }' ${de_results/DE_results/}forMapMan.txt |sed -r "s/\^\S+//g" >> ${de_results/DE_results/}forMapMan.qval0.05.txt

		[[ -z ${de_results+x} ]] || perl $wdir/process_ASE_DE_results.pl -d $de_results -t $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix -g $group_list -uniq -grp_cov 100 -min_read 50  -p1 L-LA-Purple -p2 L-US56-14 -oin L-9202 -oin S-LA-Purple -oin S-US56-14 -oin S-9202

        cd $de_fold
        for gspt in ${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix.group_uniq.minread50.grpCov100.table.*.table; do
        [[ -z ${de_results+x} ]] || perl /home/ratnesh.singh/Scripts/Perl/process_mapman_mercator.pl -gexp $de_fold/${gspt} -l 4 -i trinity -d ":" -c 1 -s 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 -o wAnnot_${gspt} -f $wdir/99_mercator_annotation_LowHighsmallerGrpFinal/LAP.US56.9202.combined.longest_isoform_as_genes.mercator.results.txt
        done

	done









	### create_gene lengths table to be used for ASE. since evey base is considered as alone, make the effective length as 1
	#$gene_length
	#[[ -e $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list ]] || cut -f1 $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix|sort|uniq > $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list
	####perl -e 'open LEN,"<","$1";while(<LEN>){my($name,$len)=split /\s+/;$seq_len{$name}=$len;} open LIST,"<","$2"; while(<LIST>){chomp;my($iso,$pos,$base)=split /:/; print "$_\t$seq_len{$iso}\n"} ' $gene_length $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list > $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list.len

	#[[ -e $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list.len ]] || while read gene ext; do printf "$gene\t1\n"; done < $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list > $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list.len


	### create GO annotation table.
	#$de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list

	#cd $de_fold && $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.EXPR.matrix --output $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.TMM.AS_DE --samples $group_list --examine_GO_enrichment --GO_annots $goannot  --gene_lengths $de_fold/${short_outfile}.bowtie2.sorted.rg.bam.freebayes.vcf.RC.OLPA.list.len -P $qval -C $log2foldc






done

cd $curdir

echo "Reached END of the Script $0 "
