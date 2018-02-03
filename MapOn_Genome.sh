#!/bin/bash
query=$1
db=$2


[[ ! -z $query ]] || { echo "$0 QUERY database " ;exit; }
[[ ! -z $db ]] || { echo "$0 query DATABASE " ;exit; }

percent_thresh=85
length_thresh=50
evalue_thresh=1E-01



mkdir "tempdb_blast2n"
ln -s $query ./tempdb_blast2n
ln -s $db ./tempdb_blast2n
cd tempdb_blast2n

#files=`ls`
#read -a file <<<$files
file1=${query##*/}
file2=${db##*/}


echo "-->makeblastdb -in $file2 -dbtype nucl"
makeblastdb -in $file2 -dbtype nucl




perl ~/Scripts/Perl/parralel_blast.pl -p blastn -n 20 -o $file1.vs.$file2.blastn.table -of '6 std qcovhsp qcovs qlen slen' -a 'culling_limit 2' -a 'num_threads 40' -a 'evalue 1e-10' -f $file1.vs.$file2.blastn -query $file1 -db $file2 -a 'xdrop_gap 500' -r


rm *.nhr *.nin *.nsq
cd ..
cp tempdb_blast2n/$file1.vs.$file2.blastn.table ./
rm -R tempdb_blast2n

perl -e 'while(<>){my$line=$_;($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,@sv)=split(/\s+/,$_);if($blast{$on}{'qcov'}<$fr){$blast{$on}{'sub'}=$tw;$blast{$on}{'qcov'}=$fr;$blast{$on}{'qlen'}=$fv;$blast{$on}{'subname'}=join " ",@sv; $blast{$on}{'all'}=$line }} foreach $query(keys %blast){print "$blast{$query}{'all'}\n"}' <$file1.vs.$file2.blastn.table>$file1.vs.$file2.blastn.TopHits.table

perl ~/Scripts/Perl/Graphics/perlgraphics_Mod.pl -i $file1.vs.$file2.blastn.TopHits.table -d subject
