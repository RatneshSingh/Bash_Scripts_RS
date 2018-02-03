### this script will extract marker names from joinmap_loc file, extracts marker sequences from tags.tsv, blasts to reference genome and cds, extracts the location of top CDS hit from GFF file. It will add the the location of top hits on the marker names and create a summary table with all the information per marker per line.

tags=$1
loc_file=$2
#### marker file with 3 colums <Markername> <distance_in_cM> <Linkage_Group>
#### eg.  56752 3.56 LG5
male_markers=$3
female_markers=$4
ref=$5

work_dir=$6

### usage
if [ -z $1 ]; then
  echo "********************"
  echo "Usage: script <catalog.tags.tsv_file> <joinMap_loc_file> <male_markers_table> <female_markers_table> <reference_species> <work_dir>"
  echo "reference options:
  Os: Oryza sativa [default]
  Sb: Sorghum bicolor
  Bd: Brachypodium distachyon
  Zm: Zea mays
  Si: Setaria italica
  At: Arabidopsis thaliana
  Pv: Panicum virgatum
  Pg: Pennisetum glaucum
  "
  echo "********************"

  exit
fi


### default work dir to current dir
if [ -z $work_dir ]; then work_dir=$PWD; fi



########

cds_database="/home/Databases/BlastDatabases"
genome_database="/home/Databases/BlastDatabases"
gff_database="/home/ratnesh.singh/Phytozome_data"

cd $work_dir
rm $ref.GenomeBlast.table

if   [ $ref -eq "Os" ];then
  echo "Osativa_204_v7.0.hardmasked.fa	Osativa_204_cds_primaryTranscriptOnly.fa  Osativa_204_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Bd" ];then
  echo "Bdistachyon_192.fa	Bdistachyon_192_cds_primaryTranscriptOnly.fa  Bdistachyon_192_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Zm" ];then
  echo "Zmays_284_AGPv3.hardmasked.fa	Zmays_284_AGPv3.transcript.fa Zmays_181_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Si" ];then
  echo "Sitalica_164_hardmasked.fa	Sitalica_164_cds_primaryTranscriptOnly.fa Sitalica_164_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Sb" ];then
  echo "Sbicolor_v2.1_255.hardmasked.fa	Sbicolor_v2.1_255_cds_primaryTranscriptOnly.fa  Sbicolor_v2.1_255_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Pv" ];then
  echo "Pvirgatum_273_v1.0.hardmasked.fa	Pvirgatum_273_v1.1.cds.fa Pvirgatum_273_v1.1.gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "At" ];then
  echo "Athaliana_167_TAIR9.fa	Athaliana_167_transcript.fa Athaliana_167_gene.gff3" >$ref.GenomeBlast.table
elif [ $ref -eq "Pg" ];then
  echo "Pglaucum.chr.fa Pglaucum.trans.cds.fa Pglaucum.trans.gff" >$ref.GenomeBlast.table
else
  echo "Using default reference:Oryza sativa"
  echo "Osativa_204_v7.0.hardmasked.fa	Osativa_204_cds_primaryTranscriptOnly.fa  Osativa_204_gene.gff3" >$ref.GenomeBlast.table
fi


if [ ! -z $male_markers ]  && [ -e "$male_markers" ];then
  #if [ ! -e "$female_markers" ]; then echo;echo "Couldnot find $female_markers";echo; exit; fi
  cur_fold="Father_RADseq_synteny"
  markers_table=$male_markers



  echo "Placing symlink in output folder $work_dir/$cur_fold"
  mkdir -p $work_dir/$cur_fold
  ln -s $tags                  $work_dir/$cur_fold/
  ln -s $loc_file              $work_dir/$cur_fold/
  ln -s $cds_database/$cds     $work_dir/$cur_fold/
  ln -s $gff_database/$gff     $work_dir/$cur_fold/
  ln -s $genome_database/$gen  $work_dir/$cur_fold/

  ln -s $markers_table         $work_dir/$cur_fold/
  ln -s $ref.GenomeBlast.table $work_dir/$cur_fold/


  cd $work_dir/$cur_fold

  while read gen cds gff ext;do
    echo "Adding reference genome and Gene model locations for $gen for Mother LG map"
    sh ~/Scripts/Bash_scripts/Add_reference_location_onmarkers.sh -c ${tags##*/} -l ${loc_file##*/} -g ${cds##*/} -G ${gff##*/} -r ${gen##*/} -f $work_dir/$cur_fold/  -F ${$markers_table##*/} -e 1e-5 -b -m; done < $ref.GenomeBlast.table

 cd $work_dir
fi


if [ ! -z $female_markers && -e "$female_markers" ];then
  #if [ ! -e "$female_markers" ]; then echo;echo "Couldnot find $female_markers";echo; exit; fi
  cur_fold="Mother_RADseq_synteny"
  markers_table=$female_markers

  mkdir -p $work_dir/$cur_fold

  echo "Placing symlink in output folder $work_dir/$cur_fold"
  ln -s $tags                  $work_dir/$cur_fold/
  ln -s $loc_file              $work_dir/$cur_fold/
  ln -s $cds_database/$cds     $work_dir/$cur_fold/
  ln -s $gff_database/$gff     $work_dir/$cur_fold/
  ln -s $genome_database/$gen  $work_dir/$cur_fold/

  ln -s $markers_table         $work_dir/$cur_fold/
  ln -s $ref.GenomeBlast.table $work_dir/$cur_fold/

  cd $work_dir/$cur_fold

  while read gen cds gff ext;do
    echo "Adding reference genome and Gene model locations for $gen for Mother LG map"
    sh ~/Scripts/Bash_scripts/Add_reference_location_onmarkers.sh -c ${tags##*/} -l ${loc_file##*/} -g ${cds##*/} -G ${gff##*/} -r ${gen##*/} -f $work_dir/$cur_fold/  -F ${$markers_table##*/} -e 1e-5 -b -m; done < $ref.GenomeBlast.table

 cd $work_dir
fi

rm GenomeBlast.table

echo;echo;
exit