#!/bin/bash

###### parsing options
function usage()
{
cat << EOF
usage: $0 options
Required options:
    -f       : Full path to fastq file containing reads.
    -b       : Barcode file containing all the barcodes. can combine barcodes of different lengths in one file.
    -o       : Output folder to save demultiplexed read files.

optional (will use defaults if not provided)
    -t       : trim all reads to length [95]
    -r       : Restriction enzymes used adapter [nsiI]
                 Currently supported enzymes include:
                 : 'apeKI', 'bamHI', 'dpnII', 'ecoRI', 'ecoT22I', 'hindIII',
                 : 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII', 'notI',
                 : 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI',
                 : 'sphI', or 'xbaI'.
    -q       : specify how quality scores are encoded [phred33],
             : 'phred33' (Illumina 1.8+, Sanger) or 'phred64' (Illumina 1.3 - 1.5).
EOF

}


while getopts “f:b:o:t:r:q:” OPTION
do
     case $OPTION in
         f)
              fastq_file=$OPTARG
              ;;
         b)

              barcode_file=$OPTARG
              ;;
         o)
              output_folder=$OPTARG
              ;;
         t)
              trim_to=$OPTARG
              ;;
         r)
              rest_enz=$OPTARG
              ;;
         q)
              phred_score=$OPTARG
              ;;

         ?)
             usage
             exit
             ;;
     esac
done







if [[ -z $fastq_file ]] || [[ -z $barcode_file ]]
then
     printf "\nERROR: fastq files of reads and barcode files are required to run the script.\n\n"
     usage
     exit 1
fi


curdir=$(pwd)
##### Default values parameters.
    : ${trim_to:=95}                      #### final trimmed length for processed reads via process_radtags.
    : ${rest_enz:="nsiI"}
    : ${phred_score:='phred33'}

#cd $output_folder

##create barcode individual files
cd $curdir
rm $barcode_file.*.base.txt
while read barcode|| [ -n "$barcode" ]; do code_len=${#barcode}; echo $barcode>>$barcode_file.$code_len.base.txt; done < $barcode_file


  base_name=`basename $fastq_file`
  ln -s  $fastq_file $output_folder/$base_name
  #base_name=`echo $fastq_file|sed 's/.fastq//'`


  for bar_files in `ls $barcode_file.*.base.txt`;do
      process_radtags -r -q -E $phred_score -t $trim_to -D --filter_illumina \
      -f $output_folder/$base_name \
      -b $bar_files \
      -o $output_folder \
      -e $rest_enz \
      #--adapter_mm 3 \
      #--adapter_1 $adapter_1
      mv $output_folder/${base_name}.discards $output_folder/$base_name
      cd $output_folder
      rename sample ${base_name%%.fastq} sample*
      ext=$(echo $bar_files|tail -c 11|head -c 6)
      rename process ${base_name%%.fastq}.${ext}.process process_radtags.log
      cd $curdir

done
  mv  $output_folder/$base_name  $output_folder/${base_name}.discards
  cd $curdir
