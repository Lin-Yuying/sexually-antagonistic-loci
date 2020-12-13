sim=1
times=2
VCF=3
OUT=4
while [ $sim -lt `$times + 1` ]
do
    echo "############ Simulation ${sim} ###############" >> pop.perm.log
    python permutation.py pop.txt male_${sim}.txt female_${sim}.txt 
    vcftools --vcf ${VCF} --weir-fst-pop male_${sim}.txt --weir-fst-pop female_${sim}.txt --out pop.snp.${sim}
    sim=`expr ${sim} + 1`
done

cat pop.snp.*.fst | grep -v CHROM > pop.snp.total.fst
cat pop.snp.*.fst | grep "CHROM"|uniq > header.txt
cat header.txt pop.snp.total.fst > pop.snp.fst
# sort SNPs according to genomic region 
python3 ./fst2csv.py pop.snp.fst pop.snp.sorted.fst 
# calculate p value for each SNP
python3 ./cal_pvalues_lin.py -in_fst pop.snp.sorted.fst -out ${OUT}
