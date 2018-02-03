gff=$1

[[ -e $gff ]] || exit


cur_dir=$PWD
[[ -d /dev/shm/split_hmmsearch ]] && rm -r /dev/shm/split_hmmsearch
mkdir -p /dev/shm/split_hmmsearch
cp $gff /dev/shm/split_hmmsearch/
cp /home/Databases/pfam/Pfam-A.hm* /dev/shm/split_hmmsearch/



if [[ ! -e ${gff%.*}.prot.fasta ]]; then
    /home/Databases/softwares/JAMg/3rd_party/evidencemodeler/EvmUtils/gff3_file_to_proteins.pl $gff $genome prot > ${gff%.*}.prot.fasta	
    
fi

prot=${gff%.*}.prot.fasta

cp $prot /dev/shm/split_hmmsearch/
cd /dev/shm/split_hmmsearch/

### convert files names to basenames

gff=$(basename $gff)
prot=$(basename $prot)
####



echo "Splitting $prot in chunks of 500 sequences"
splitfasta.pl -i $prot -depth 500 -suffix split -dir .


echo "Running hmmsearch against PfamA "
for file in *.split; do echo hmmsearch  -o ${file}.hmmscan.table --tblout ${file}.hmmscan.tblout --acc --noali -E 0.1 --cpu 3 Pfam-A.hmm $file ;done |parallel -j 63

echo "Combining results from all the hmmsearch runs"
cat *.hmmscan.table > ${prot%.*}.hmmscan.table2

cat *.hmmscan.tblout > ${prot%.*}.hmmscan.tblout2

echo "Removing split files and their output"
rm *.split* 

echo "Creating list of Sequences showing Pfam domains"
grep -v "^#" ${prot%.*}.hmmscan.tblout2|cut -f1 -d " "|sort|uniq > ${prot%.*}.hmmscan.tblout2.list


perl ~/Scripts/Perl/search_seq-V3.4.pl -s ${prot} -l ${prot%.*}.hmmscan.tblout2.list -o ${prot%.*}.hmmscan.tblout2.fasta -m exact -d " " -c 1

echo "Removing Repeat containing sequences"
echo "removing the Pfam related to Repeats from hmmsearch results"
egrep "PF00075|PF00077|PF00078|PF00098|PF00424|PF00429|PF00469|PF00516|PF00517|PF00522|PF00539|PF00540|PF00552|PF00558|PF00559|PF00589|PF00607|PF00665|PF00692|PF00872|PF00906|PF00971|PF00979|PF01021|PF01045|PF01054|PF01140|PF01141|PF01359|PF01385|PF01498|PF01526|PF01527|PF01548|PF01609|PF01610|PF01695|PF01710|PF01797|PF02022|PF02093|PF02228|PF02281|PF02316|PF02337|PF02371|PF02411|PF02720|PF02813|PF02892|PF02914|PF02920|PF02959|PF02992|PF02994|PF02998|PF03004|PF03017|PF03050|PF03056|PF03078|PF03108|PF03184|PF03221|PF03274|PF03276|PF03400|PF03408|PF03539|PF03708|PF03716|PF03732|PF03811|PF04094|PF04160|PF04195|PF04218|PF04236|PF04582|PF04693|PF04740|PF04754|PF04827|PF04937|PF04986|PF05052|PF05344|PF05380|PF05399|PF05457|PF05485|PF05598|PF05599|PF05621|PF05699|PF05717|PF05754|PF05840|PF05851|PF05858|PF05928|PF06527|PF06815|PF06817|PF07253|PF07282|PF07567|PF07572|PF07592|PF07727|PF07999|PF08284|PF08333|PF08483|PF08705|PF08721|PF08722|PF08723|PF09035|PF09039|PF09077|PF09293|PF09299|PF09322" ${prot%.*}.hmmscan.tblout2|cut -f1 -d " "|sort|uniq > ${prot%.*}.hmmscan.tblout2.WidRepeatPfam.list


perl ~/Scripts/Perl/Remove_selectedSequences-V1.3.pl -s ${prot%.*}.hmmscan.tblout2.fasta -l ${prot%.*}.hmmscan.tblout2.WidRepeatPfam.list -m exact -d " " -c 1
 
 
echo "Copying resulting fasta and gff files back to $cur_dir"
cp ${prot%.*}.hmmscan.tblout2.fasta.include.fasta $cur_dir/${prot%.*}.NoPfamARepeat.fasta 
cp ${prot%.*}.hmmscan.tblout2.fasta.removed.fasta $cur_dir/${prot%.*}.ContainsPfamARepeat.fasta
cp ${prot%.*}.hmmscan.tblout2 $cur_dir/${prot%.*}.hmmscan.tblout2

echo "Filtering GFF3 file for sequences containing Pfam Domains but not Repeats"
grep ">"   ${prot%.*}.hmmscan.tblout2.fasta.include.fasta |sed 's|>||g' > ${prot%.*}.NoPfamARepeat.list
grep -F  -w -f ${prot%.*}.NoPfamARepeat.list $gff > $cur_dir/${gff%.*}.NoPfamARepeat.gff3


cd $cur_dir/




echo "Running busco for completenes test on ${prot%.*}.NoPfamARepeat.fasta"
busco -f -c 45 -l /mnt/nfs1/usr/local/busco/datasets/embryophyta_odb9  -i ${prot%.*}.NoPfamARepeat.fasta -o run_${prot%.*}.NoPfamARepeat.busco_maize_out -m prot  &



echo "Finished running $0"

echo
