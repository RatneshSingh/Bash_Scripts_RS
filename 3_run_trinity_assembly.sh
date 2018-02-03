### Remember to run S-7066 seperately as it was giving problem due to missing R2 read files.
### read file names and collect names of sample. Run trinity for each sample (leaf + stem combined)
## Assemblies done:
## _mkc2_mpisp90_mdsp8_migsp10
## _mkc2_mpisp95_mdsp5_migsp15
## _mkc2_mpisp95_mdsp8_migsp10
## _mkc2_mpisp95_mdsp8_migsp10
## _mkc2_mpisp97_mdsp10_migsp15
## _mkc4_mpisp90_mdsp8_migsp10_kmer31
## 

mkc=2     # --min_kmer_cov      
mpisp=97  # --min_per_id_same_path  <int> 
mdsp=8    # --max_diffs_same_path <int>
migsp=10  # --max_internal_gap_same_path <int>
kmer=31   # --KMER_SIZE 

trinity="/usr/local/trinity/Trinity"

for name in $(ls /home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/1_raw_reads/|grep "fastq"|cut -f1 -d "_"|cut -f2,3 -d "-"|sort|uniq|egrep "(^U|^92|^L)"); do
  
  	out_fold="3_trinity_assemblies/trinity.${name}_mkc${mkc}_mpisp${mpisp}_mdsp${mdsp}_migsp${migsp}_kmer${kmer}"
  	#mkdir -vp $out_fold

  	#ln -s $PWD/2_trimmed_reads/$name/*.fq $PWD/2_trimmed_reads/$name/*.gz $PWD/2_trimmed_reads/$name/*.readcount    $PWD/2_trimmed_reads/$name/both*   $PWD/$out_fold/
	rm -r $PWD/$out_fold
	mkdir -p $out_fold
	for file in `ls $PWD/2_trimmed_reads/$name/*.{fa,fq,gz,ok,count}`; do echo ln -fs $file $PWD/$out_fold/ ; done

  	## Read all the file for sample $name
  	R1=$(ls /home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/1_raw_reads/*R1.fastq.gz|grep $name);
  	R1_seqs=${R1[@]};

  	left=${R1_seqs// /,}
 	 right=${left//R1.fastq/R2.fastq}
 	  #echo "--left $left --right $right"

	  printf  "$trinity --KMER_SIZE $kmer --full_cleanup  --seqType fq --max_memory 300G --SS_lib_type RF --CPU 30 --trimmomatic --output $out_fold  --min_kmer_cov  $mkc --min_per_id_same_path $mpisp --max_diffs_same_path $mdsp --max_internal_gap_same_path $migsp  --verbose  --left ${left//$'\n'/} --right ${right//$'\n'/} > ${out_fold##*/}.log 2>&1  \n"



done > Trinity_commands_to_run.list

dos2unix Trinity_commands_to_run.list


### Run all the commands using parallel command.
### -j 5 : run 5 commands at a time.
parallel  -j 10 < Trinity_commands_to_run.list 

