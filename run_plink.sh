## usage: script plink file address    out folder address  folder_name for plink file (.ped, .map)
curdir=$(pwd)


## ped file full path
plink_file=$1

plink_name=${plink_file##*/}
## out folder full path
folder=$2;



#rm -R $folder/plink
mkdir -p $folder/plink

#plink_file=$(ls|grep ".ped$")

cd $folder/plink
## convert plink result to binary file
plink --allow-no-sex --file ${plink_file%%.ped} --make-bed --out $folder/plink/${plink_name%%.ped}
cd $folder/plink
## update names as new names without species anems or extra extensions
# while read fid iid ext; do niid=${iid%%[_|Q]*};echo $fid $iid $fid ${niid#*.}; done< ../$plink_file|grep -v '#'  >${plink_file%%.ped}name.update

plink --allow-no-sex  --file ${plink_file%%.ped} --out ${plink_name%%.ped} --update-ids ${plink_name%%.ped}name.update  --make-bed

plink --allow-no-sex  --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --recodeA

#####generate some simple summary statistics on rates of missing data in the file, using the --missing option
## genotyping rate
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --missing

#####The following command generates a file called freq_stat.frq which contains the minor allele frequency and allele codes for each SNP.
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --freq

## number monomorphic sites.
## filters:  --mind individual missing rate; --geno =genotyp missing rate; --maf= SNP allele frequency; maximum-allele-frequency
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped}  --mind 1 --geno 1 --maf 0 --max-maf 0


## non random genotyping failure
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --test-missing

## create binay file from filtered data
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped}.filtered --make-bed --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-3

## Most highly associated SNP. Out: * .assoc (sorted by physical position) and *.assoc.adjusted (sorted by p-value)
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped}.filtered  --assoc --adjust



#### perform a cluster analysis that pairs up individuals on the basis of genetic identity.which requests IBS clustering (--cluster) but with the constraints
#### that each cluster has no more than two individuals (--mc 2) and that any pair of individuals who have a significance value of less than 0.05 for the test
#### of whether or not the two individuals belong to the same population based on the available SNP data are not merged.
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --cluster


####This is a simple significance test for whether two individuals belong to the same random-mating population.
#### To only merge clusters that do not contain individuals differing at a certain p-value (smaller the p-value, smaller the group number) -ppc:
plink --noweb --allow-no-sex --bfile ${plink_name%%.ped} --cluster --ppc 0.0001 -mc 2 --out ${folder}.plink.ppc0.0001


#### Given that pairwise IBS distances between all individuals have been calculated, we can asked whether or not there are
#### group differences in this metric, with respect to a binary phenotype.
plink --noweb --allow-no-sex --bfile ${plink_name%%.ped} --ibs-test --out ${plink_name%%.ped}


#### calculate genome-wide IBD given IBS information, as long as a large number of SNPs are
#### available (probably 1000 independent SNPs at a bare minimum; ideally 100K or more).
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --genome


#### Given a large number of SNPs, in a homogeneous sample, it is possible to calculate inbreeding
#### coefficients (i.e. based on the observed versus expected number of homozygous genotypes).
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --het


##### For the N individuals in a sample, to create a N x N matrix of genome-wide average IBS pairwise identities
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --cluster --matrix


#### For the N individuals in a sample, to create a N x N matrix of genome-wide average IBS pairwise distance identities
plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped} --cluster --distance-matrix



## Multidimensional scaling plots
plink --file mydata --read-genome ${plink_name%%.ped}.genome --cluster --mds-plot 4 




## write covariate file
#plink --allow-no-sex --bfile ${plink_name%%.ped} --out ${plink_name%%.ped}.filtered --write-covar myfile.txt




####===============================================================================

cd $curdir
