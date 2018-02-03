seqfile=$1

cur_dir=$PWD
mkdir -p ${seqfile%.*}_ACT

cd ${seqfile%.*}_ACT

perl -e '$/="\n>"; while(<>){my($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; my$shead=$header;$shead=~s/\s+.*$//g;open OUT, ">$shead.fa";  print OUT ">$header\n$sequence\n"; close OUT;} ' < $cur_dir/$seqfile


shopt -s nullglob
array=(*.fa)
larray=$(( ${#array[@]} - 1 ))
i=0

while [[ i -lt ${#array[@]} ]] ; do
  j=$((i+1))
  while [[ j -lt ${#array[@]} ]] ; do
    echo "Running blast between ${array[$i]} ${array[$j]}"
    sh /home/ratnesh.singh/Scripts/Bash_scripts/blast2dna.sh ${array[$i]} ${array[$j]}
    echo "done"
    j=$((j+1))
  done
  i=$((i+1))
done

cd $cur_dir