upfile=$1;
downfile=$2
#seq_fasta=$1

for seq_fasta in $upfile $downfile; do
  echo "*** Running the KAAS annotation on the web server for $seq_fasta"
  #### run KAAS to annotate sequences and map to pathways. Save the jobID
  sed -i 's/[ \t]*//g' ${seq_fasta/fasta/}jobid
  jobid=$(cat ${seq_fasta/fasta/}jobid|sed 's/[ \t]*//g')

  if [[ $jobid == "" ]]; then
      echo "*** Running the KAAS annotation on the web server for $seq_fasta"
      curl -X POST  --form continue=1  --form prog=BLAST  --form uptype=q_text  --form text=  --form peptide2=n  --form "file=@${seq_fasta}"  --form qname=${seq_fasta/.fasta/}  --form mail=ratnesh.singf@ag.tamu.edu  --form dbmode=manual  --form "org_list=ath, aly, tcc, gmx, fve, csv, vvi, sly, osa, olu, ota, mis, cme, gsl"  --form way=s  --form mode=compute http://www.genome.jp/kaas-bin/kaas_main > ${seq_fasta/fasta/}KAAS_jobsubmission.log

      jobid=$(grep  "job ID" ${seq_fasta/fasta/}KAAS_jobsubmission.log|cut -f 2 -d ":"|cut -f 1 -d "<"|sed 's/^[ \t]*//;s/[ \t]*$//')

      echo $jobid > ${seq_fasta/fasta/}jobid
      echo "The jobid of submitted jb is : $jobid"

  else
      echo "jobid:$jobid exists. Not running new task. Just going to retrieve the file from server."
  fi


  #jobid=1437078203
  ##### retrieve the result as KO file from KAAS server
    if [[  -s ${seq_fasta/fasta/KAAS}.ko ]]; then

      echo "${seq_fasta/fasta/KAAS}.ko exists. Please delete it or rename to fetch again"
      break
    else

      [[ $jobid == "" ]] && echo "Jobid is empty: existing" && exit
      while true; do
        sleep 60
        wget --spider "http://www.genome.jp/kaas-bin/kaas_main?mode=map&id=${jobid}&mail=ratnesh.singf@ag.tamu.edu"
        check=$(wget --spider http://www.genome.jp/tools/kaas/files/dl/$jobid/query.ko|grep broken)

       [[ "$check"  =~  "" ]] &&  wget -O ${seq_fasta/fasta/KAAS}.ko http://www.genome.jp/tools/kaas/files/dl/$jobid/query.ko && break
      done
    fi

done



#### filter the KO number and assign red and green color for mapping through kegg mapper.
if [[ -s  ${upfile/fasta/KAAS}.ko ]] && [[ -s  ${downfile/fasta/KAAS}.ko ]]; then
grep -P "\s+K" ${upfile/fasta/KAAS}.ko|cut -f 2|sed 's/[ \t]*$/ red/' > ${upfile/up/UpDown}.color.list
grep -P "\s+K" ${downfile/fasta/KAAS}.ko|cut -f 2|sed 's/[ \t]*$/ green/' >> ${upfile/up/UpDown}.color.list

file=${upfile/up/UpDown}.color.list
fi

[[ $file == "" ]] && echo "color list file:$file seems as empty name" && exit
############# upload color list data and get mappped page back
[[ -f ${file/list/htm} ]]  || curl -X POST --form "color_list=@${file}" --form org=ko --form s_sample="" --form unclassified="" --form default=pink --form target=alias --form reference=white --form warning=yes  http://www.genome.jp/kegg-bin/color_pathway_object > ${file/list/htm}

### wait 5 seconds
sleep 5
##### process the received result from KEGG mapper to retrieve mapped pathways
##### and retract marked png.
folder=${file/.list/_mapped}
[[ -d $folder ]] && exit
mkdir -p $folder
cd $folder
grep  'args" target=' ../${file/list/htm} > ${file/list/htm}.pathway.links
 while read link ext; do
   alink=${ext#*\"};
   plink=${alink%%\"*}
   wget  http://www.genome.jp${plink};
done < ${file/list/htm}.pathway.links


 for png in *.args; do
   link=$(grep .png\"\ name=\"pathwayimage\" $png);
   alink=${link#*\"};
   nlink=${alink%%\"*};
   [[ -f ${file%.*}.${nlink##*/} ]] || wget -O ${file%.*}.${nlink##*/}  http://www.genome.jp${alink%%\"*};

  done
