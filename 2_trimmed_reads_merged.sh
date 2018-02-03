cur_dir=$PWD
readdir="2_trimmed_reads"
raw_fold="1_raw_reads"
destdir="2_trimmed_reads_merged"
for fold in $cur_dir/$readdir/*;do 

	#echo $fold
	R1files=$(ls  $fold/*R1*P.qtrim.gz|cut -d "_" -f 6,7|sort|uniq)
	dfold=${fold##*/}
	destfold="$cur_dir/$destdir/$dfold"

	for file in $R1files; do
	   if [ -f $raw_fold/${file##*/}*R1.fastq.gz ] && [ -f $raw_fold/${file##*/}*R2.fastq.gz ]; then 
              if [ -f $destfold/${file##*/}.R1n2.qtrim.merged.assembled.fastq ]; then
	        echo "output file $destfold/${file##*/}.R1n2.qtrim.merged.assembled.fastq  from previous run exists. Please delete it to rerun for this file set"
              else	
		mkdir -p $destfold
                file_f=$(ls $raw_fold/${file##*/}*R1.fastq.gz)
		file_r=$(ls $raw_fold/${file##*/}*R2.fastq.gz)
 		pearRM  -j 50 -f $file_f -r $file_r -o $destfold/${file##*/}.R1n2.qtrim.merged
              fi
	   else
 		echo "Both pairs do not exist"
           fi
        done 

done 
#> 2_trimmed_reads_merged.command.list
