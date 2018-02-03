assembly=$1

[[ -f $assembly ]] || (echo "Require assembly name to test"; exit)

asname=${assembly##*/}

#for name in $(ls /home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/1_raw_reads/|grep "fastq"|cut -f1 -d "_"|cut -f2,3 -d "-"|sort|uniq|egrep "(^U|^92|^L)"); do
for name in $(echo $asname|cut -f 1 -d "_"|cut -f 2 -d '.'); do  
  	out_fold="/home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/6_evaluate_assemblies/detonate.${name}"

	mkdir -p $out_fold
slen=$(head -n 2 /home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/2_trimmed_reads/$name/both.fa|tail -n1 |wc -m)

rsem-eval-calculate-score -p 20 \
              --transcript-length-parameters /mnt/nfs1/ratnesh.singh/.usr/local/detonate-1.10/rsem-eval/true_transcript_length_distribution/Sbicolor.txt \
              /home/ratnesh.singh/Sugarcane/RNA_Seq_Hawaii_populations/2_trimmed_reads/$name/both.fa \
              $assembly  \
              $out_fold/${asname/fasta/}rsem_eval  \
              $slen  

done 

