function usage() {
 echo "$0 ...options
 Options:
    -g reference genome in fasta format
    -n gene name to draw splice graph for.
    -P use mRNA name from gff .
    -o outfile name to save results [$outfile ]
    -c classifier file [will be created from gff]
    -h print this usage and exit.

  SAM alignment file. Provide atleast one of these to show Junctions and Read coverage tracks.
    -s unsorted sam file containing alignments. will be filtered using sam_filter.py
    -S sorted sam file containing alignments. will be filtered using sam_filter.py
    -f filtered sorted sam file.
    
  EST related files. Provide atleast one to show EST pannel in figure.
    -e est fasta file. Will be used in blat alignment to generate PSL and est_gff.
    -p blat alignment of ests in psl format. Will be used to generate est_gff file using ests_to_splicegraph.py.
    -G EST gff from ests_to_splicegraph.py tool. Will be used directly in drawing.
 
  Gene structure related. Provide atleast one of these.
    -r mRNA in fasta format
    -m gene model in gff format. mRNA file will be used to create one if not provided.
    -P use gene's name in gff3 file to draw figure for.

  Splicegrapher options
    -J Minimum depth required for junction evidence [1]
    -M Minimum anchor size required for junction evidence [5]
    -T Minimum depth threshold for identifying clusters [1]

    "

    exit;
}


#function main(){

#local  OPTIND OPTARG sam  ssam fsam est psl genome est_gff mRNA gene_model gene_name classifier J M T outfile gm_as_name



while getopts ":s:S:f:e:p:g:G:r:m:n:c:J:M:T:o:hP" opt; do
    case $opt in
       s )  sam=$OPTARG ;;
       S )  ssam=$OPTARG ;;
       f )  fsam=$OPTARG ;;
       e )  est=$OPTARG ;;
       p )  psl=$OPTARG ;;
       g )  genome=$OPTARG ;;
       G )  est_gff=$OPTARG ;;
       r )  mRNA=$OPTARG ;;
       m )  gene_model=$OPTARG ;;
       n )  gene_name=$OPTARG ;;   
       c )  classifier=$OPTARG ;;
       J )  J=$OPTARG ;;
       M )  M=$OPTARG ;;
       T )  T=$OPTARG ;;
       o )  outfile=$OPTARG ;;
       P )  gm_as_name=1;;
       h )  usage ;;
       \?)  echo "Unknown argument $OPTARG"; usage ;;
   esac
done
###########################################################################################################
#### create classifier from gff and genome if not provided by the user
if [[ -z "$classifier" && ! -z "$genome"  && ! -z "$gene_model" ]]; then
    echo "No classifier is provided. Program will try to create one"
	if [[ -f $genome  && -f  $gene_model ]]; then
		build_classifiers.py -d gt,gc,at -a ag,ac -f $genome -m $gene_model
		classifier="classifiers.zip"
	else
        echo "Program cannot create classifier since genome and gff file is not provided. Will use default arabidopsis's classifier"
    fi
else
	echo "Will use classifier: $classifier"
	
fi

: ${classifier:=/home/Databases/softwares/SpliceGrapher-0.2.4/classifiers/Arabidopsis_thaliana.zip}

#############################################################################################################

: ${J:=1}    ## Minimum depth required for junction evidence
: ${M:=5}   ## Minimum anchor size required for junction evidence
: ${T:=1}    ## Minimum depth threshold for identifying clusters
: ${classifier:=/home/Databases/softwares/SpliceGrapher-0.2.4/classifiers/Arabidopsis_thaliana.zip}
[[ -z $gene_model ]] || gff_name=$(basename $gene_model)

if [[ -z "$sam" && -z "$fsam" && -z "$ssam" ]]; then 
    echo "SAM alignments are not provided. Junctions and Reads will not be included in final figure."; 
    
fi 


###=========================================================================================================
###== gene model pannel
if [[ -z "$gene_model" && -z "$mRNA" ]]; then
    echo "Gene model or mRNA was not provided. please provide atleast one of two."
    usage
else
   if [[ -z  "$gene_model" ]]; then
     gmap -g $genome -f gff3_gene $mRNA > ${mRNA%.*}_gmap.gff3
     
     gene_model=${mRNA%.*}_gmap.gff3
     sed -i 's|.path1||g' $gene_model
      
   fi
    ## use gene name from gene model file if not provided by the user.
    if [[ $gm_as_name == 1 ]]; then
        mname=$(grep -v "^#" ${gene_model}|grep -w  gene |head -n 1|cut -f 9|cut -f 1 -d ';'|sed 's|ID=||g')
        gene_name=$mname
    fi

gmod_cfg="[GeneModelGraph]
plot_type     = gene
gene_name     = $gene_name
relative_size = 10.0
source_file   = $gene_model
file_format   = gene_model
title_string  = Gene Model for %gene
"
fi

### create outfile name if not provided.
: ${outfile:=${gene_name}_splicegrapher.pdf}

if [[ ! -z "$gene_name" ]]; then
    gene_model_to_splicegraph.py -m $gene_model -g $gene_name -a -o ${gene_name}_graph.gff
fi


####==========================================================================================================
#######== filter sam alignments to remove noise and nonsplicing elements.

if [[ -z "$fsam" ]] ; then
  if [[ -z "$ssam" ]]; then
    if [[ ! -z "$sam" ]]; then
       samtools sort -@ 6 -o ${sam%.*}sorted.sam $sam
       ssam=${sam%.*}sorted.sam
    fi
  fi
  
  if [[ ! -z "$ssam" ]]; then
   	## filter sorted samfile
   	sam_filter.py $ssam $classifier -f $genome -m $gene_model -o ${ssam%.*}.filtered.sam
   	fsam=${ssam%.*}.filtered.sam
   	
  fi
fi

fsam_name=$(basename $fsam)




### convert sam files to depth file
if [[ -f $fsam && ! -f ${fsam%.*}.depths ]]; then
	sam_to_depths.py $fsam -o ${fsam%.*}.depths 
	fsam_depth=${fsam%.*}.depths
fi





if [[ ! -z  "$fsam" ]]; then
sam_cfg="[Junctions]
plot_type     = junctions
labels        = True
relative_size = 5.0
source_file   = $fsam
title_string  = Splice Junctions


[Reads]
plot_type     = read_depth
relative_size = 8.0
source_file   = $fsam
title_string  = Read Coverage
"

fsam_cmd="-d $fsam"


fi



###==========================================================================================================
###== EST pannel

if [[ -z "$psl" && -z "$est" && -z "$est_gff" ]]; then
    echo "EST or PSL file is not provided. Will not add ESTS pannel to figure."
else
    if [[ -z "$est_gff" ]]; then
        if [[ -z "$psl" ]]; then
               blat $genome $est ${est%.*}.psl
               psl=${est%.*}.psl
               sed -i '/^\s*match/d;/^-/d;/^psLayout/d;/^\s*$/d' $psl
               ests_to_splicegraph.py $psl  -i 10 -g $gene_name  -m $gene_model
               est_gff=${gene_name}_ests.gff
           
        else
               sed -i '/^\s*match/d;/^-/d;/^psLayout/d;/^\s*$/d' $psl
               ests_to_splicegraph.py $psl  -i 10 -g $gene_name  -m $gene_model
               est_gff=${gene_name}_ests.gff
        fi
    fi 

[[ -e ${gene_name^^}_ests.gff ]] && mv ${gene_name^^}_ests.gff ${gene_name}_ests.gff
est_cmd="-s $est_gff"
est_cfg="[ESTs]
plot_type     = splice_graph
relative_size = 10.0
source_file   = $est_gff
title_string  = EST evidences
"    
   

fi











####============================================================================================================
###== Generating predictions ==
if [[ ! -z "$gene_name" ]]; then
predict_splicegraph.py ${gene_name}_graph.gff $fsam_cmd  $est_cmd  -o ${gene_name}_predicted.gff -J $J -M $M -T $T
splice_cfg="[SpliceGrapher]
plot_type     = splice_graph
relative_size = 10.0
source_file   = ${gene_name}_predicted.gff
title_string  = Splice Prediction
"
fi
### predict splice graph for whole sam slignment for statistics purpose
#
#if [[ -f ${fsam%.*}.depths ]];then
#	predict_graphs.py ${fsam%.*}.depths  -m $gene_model -d ${fsam_name%.*}_predicted_graphs  -j $J -a $M -D $T 
#else
#	predict_graphs.py $fsam -m $gene_model -d ${fsam_name%.*}_predicted_graphs  -j $J -a $M -D $T
#fi

##create splice graph from SAM alignment gene model
predict_graphs.py $fsam -m $gene_model -d ${fsam_name%.*}_predicted_graphs  -j $J -a $M -D $T
if [[ -d ${fsam_name%.*}_predicted_graphs ]];then
    #predict_graphs.py $fsam  -d ${fsam_name%.*}_predicted_graphs  -j $J -a $M -D $T
    find ${fsam_name%.*}_predicted_graphs -name "*.gff" > ${fsam_name%.*}_predicted_graphs.lis
    splicegraph_statistics.py ${fsam_name%.*}_predicted_graphs.lis -o ${fsam_name%.*}_predicted_graphs.stat
    genewise_statistics.py ${fsam_name%.*}_predicted_graphs.lis > ${fsam_name%.*}_predicted_graphs.AS_summary.out
    genewise_statistics.py ${fsam_name%.*}_predicted_graphs.lis -C > ${fsam_name%.*}_predicted_graphs.AS_summary.csv
fi

## create splice graph from gene model
gene_model_to_splicegraph.py -m $gene_model -a -S -d ${gff_name%.*}_predicted_graphs
if [[ -d ${gff_name%.*}_predicted_graphs ]];then
    find ${gff_name%.*}_predicted_graphs -name "*.gff" > ${gff_name%.*}_predicted_graphs.lis
    splicegraph_statistics.py ${gff_name%.*}_predicted_graphs.lis -o ${gff_name%.*}_predicted_graphs.stat
    genewise_statistics.py ${gff_name%.*}_predicted_graphs.lis > ${gff_name%.*}_predicted_graphs.AS_summary.out
    genewise_statistics.py ${gff_name%.*}_predicted_graphs.lis -C > ${gff_name%.*}_predicted_graphs.AS_summary.csv
fi
## create stat for splice graphs in gene model and sam based predcition







if [[ ! -z "$gene_name" ]]; then

#########################################################
##### create final gene graph input cfg to draw figure
echo "[SinglePlotConfig]
legend        = True
output_file   = $outfile
width         = 9.0
height        = 11.0

$gmod_cfg

$est_cfg

$splice_cfg

$sam_cfg
" > ${gene_name}.cfg


###== Generating plots
plotter.py ${gene_name}.cfg


fi






#}


#main



