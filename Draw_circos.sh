wor_dir=$PWD
output_folder=$3
if [[ -z "$1" ]]; then
  echo "Usage script <link_file> <Conf_file> <output_folder>";
  echo "link file contains columns
<LikageGroup> <mark_start> <mark_end> <Chr_Name> <Chr_start> <Chr_end>

===========eg Example link file =======
LG1     58 59 Os3 6516807 6526807
LG1     17 18 Os3 29317674 29327674
LG1     35 36 Os7 28450085 28460085
LG1     58 59 Os3 6735251 6745251
=======================================

To run circos, do following:
1: Make sure the conf file has the name of the bundle link file in the links section.


exit;
"
fi

link_file=$1
conf_file=$2
if [[ -z "$out_folder" ]]; then output_folder="circos_figures"; fi
mkdir -p $output_folder

## save the files used in circos output folder
cp $link_file $output_folder/
cp $conf_file $output_folder/

link_name=${link_file##*/}
conf_name=${conf_file##*/}

#cd $wor_dir/circos_figures
/home/ratnesh.singh/softwares/circos-tools-0.19/tools/bundlelinks/bin/bundlelinks -links $link_file -max_gap 10000 >$wor_dir/${link_name/txt/bundle.txt}
cp $wor_dir/${link_name/txt/bundle.txt} $output_folder/
echo;echo "Make sure the link file: $wor_dir/${link_name/txt/bundle.txt}  name is added in $conf_file";echo


sed -r "s/file\s*=\s*\S+/file=${link_name/txt/bundle.txt}/" $conf_file > ${conf_file/conf/mod.conf}
cp ${conf_file/conf/mod.conf} $output_folder/

### keep the full path of circos as all the required files are in the folder.
/mnt/yulab_16TB/Databases/softwares/circos-0.67-1/bin/circos \
      -conf   $output_folder/${conf_file/conf/mod.conf} \
      -outputdir $output_folder \
      -outputfile ${link_name%.*}.${conf_name/conf/png}
