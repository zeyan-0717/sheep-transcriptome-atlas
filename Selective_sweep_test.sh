#### Selective sweep test based on Fst ####
for i in {1..27}
do
bcftools view <SET_YOUR_PATH_TO_VCF/>RD.chr${i}.snp.vcf.gz -c 1 -S <SET_YOUR_PATH>/sample.list --threads 4 -Oz > <SET_YOUR_PATH>/TIB.chr${i}.snp.vcf.gz
done 

for i in {1..27}
do
bcftools index -t <SET_YOUR_PATH>/TIB.chr${i}.snp.vcf.gz
done 


for i in {1..27}
do 
<SET_YOUR_PATH_TO_SOFTWARE>/gatk-4.2.3.0/gatk SelectVariants -restrict-alleles-to BIALLELIC -V <SET_YOUR_PATH>/TIB.chr${i}.snp.vcf.gz -O <SET_YOUR_PATH>/f1.TIB.chr${i}.snp.vcf.gz
done 


# filtering out MAF (Minor allele frequency) < 0.01 and Missing rate>0.05

for i in {1..27}
do 
<SET_YOUR_PATH_TO_SOFTWARE>/vcftools/bin/vcftools --gzvcf <SET_YOUR_PATH>/f1.TIB.chr${i}.snp.vcf.gz --maf 0.01 --max-missing 0.95 --recode --out <SET_YOUR_PATH>/f2.TIB.chr${i}.snp.vcf.gz
done


for i in {1..27}
do 
<SET_YOUR_PATH_TO_SOFTWARE>/vcftools/bin/vcftools --vcf <SET_YOUR_PATH>/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf --weir-fst-pop <SET_YOUR_PATH>/HA_sheep.list --weir-fst-pop <SET_YOUR_PATH>/LA_sheep.list  --out HA_vs_LA.chr${i} --fst-window-size 10000 --fst-window-step 10000
done

for i in {1..27}
do 
<SET_YOUR_PATH>/bin/htslib-1.12/bin/bgzip <SET_YOUR_PATH>/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf -@ 6
done 

for i in {1..27}
do 
<SET_YOUR_PATH_TO>/bin/htslib-1.12/bin/tabix  <SET_YOUR_PATH>/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf.gz 
done

for i in {1..27}
do
<SET_YOUR_PATH>/bin/Python3.8/bin/python3.8 <SET_YOUR_PATH_TO_SOFTWARE>/genomics_general-master/VCF_processing/parseVCFs.py -i <SET_YOUR_PATH>/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf.gz --threads 6 | bgzip --threads 6 > <SET_YOUR_PATH>/f2.TIB.chr${i}.geno.gz
done

for i in {1..27}
do
<SET_YOUR_PATH_TO>/bin/Python3.8/bin/python3.8 <SET_YOUR_PATH_TO_SOFTWARE>/genomics_general-master/popgenWindows.py --windSize 10000 --stepSize 10000 --popsFile <SET_YOUR_PATH>/popsFile  --outFile <SET_YOUR_PATH>/FST.chr${i}
done

for i in {1..27}
do 
bcftools view -e 'F_MISSING > 0.05 & MAF <0.01' --threads 6 <SET_YOUR_PATH>/f1.TIB.chr${i}.snp.vcf.gz -Oz -o <SET_YOUR_PATH>/f3.TIB.chr${i}.snp.vcf.gz
done 

for i in {1..27}
do
<SET_YOUR_PATH_TO_SOFTWARE>/bcftools stats f1.TIB.chr${i}.snp.vcf.gz > f1.TIB.chr${i}.stats
done

for i in {1..27}
do
cat -n f1.TIB.chr${i}.stats | grep 'number of SNPs: ' >> 1_27_stats
done

for i in {1..27}
do 
bcftools index --threads 6 <SET_YOUR_PATH>/f3.TIB.chr${i}.snp.vcf.gz
<SET_YOUR_PATH>/bin/htslib-1.12/bin/tabix  <SET_YOUR_PATH>/f3.TIB.chr${i}.snp.vcf.gz
done 


for i in {1..27}
do
source ~/.bashrc
<SET_YOUR_PATH>/bin/Python3.8/bin/python3.8 <SET_YOUR_PATH_TO_SOFTWARE>/genomics_general-master/VCF_processing/parseVCFs.py -i <SET_YOUR_PATH>/f3.TIB.chr${i}.snp.vcf.gz --threads 6 --include ${i} | bgzip --threads 6 > <SET_YOUR_PATH>/f2.TIB.chr${i}.geno.gz
done

for i in {1..27}
do
source ~/.bashrc
<SET_YOUR_PATH>/bin/Python3.8/bin/python3.8 <SET_YOUR_PATH_TO_SOFTWARE>/genomics_general-master/popgenWindows.py -g <SET_YOUR_PATH>/f2.TIB.chr${i}.geno.gz -f phased -T 5 --windSize 10000 --stepSize 10000 -p HA GNM327,GNM330,GNM341,CJL285,CJL289,CJL300,SAB311,SAB319,SAB322,ZJZ225,ZJZ226,ZJZ242,ZLX154,ZLX160,ZLX161,CHA02,CHA05 -p LA HU62,HU65,HU69,HU63,SXW30,SXW19,SXW29,SXW16,WDS123,WDS14,WDS109,WDS112,SSS44,SSS45,SSS46,SSS32,HDW47,HDW51,HDW52,HDW48  --popsFile <SET_YOUR_PATH>/popsFile  --outFile <SET_YOUR_PATH>/FST.chr${i}
done

for i in {1..27}
do
zcat f2.TIB.chr${i}.geno.gz | wc -l
done


#### Sliding window approach (10 kb sliding windows with 10 kb steps)
cat FST.chr1	FST.chr2	FST.chr3	FST.chr4	FST.chr5	FST.chr6	FST.chr7	FST.chr8	FST.chr9	FST.chr10	FST.chr11	FST.chr12	FST.chr13	FST.chr14	FST.chr15	FST.chr16	FST.chr17	FST.chr18	FST.chr19	FST.chr20	FST.chr21	FST.chr22	FST.chr23	FST.chr24	FST.chr25	FST.chr26	FST.chr27 > FST.popgenWindows

sed '^/scaffold/d' FST.popgenWindows > tmp1.FST.popgenWindows
sed '1i\scaffold,start,end,mid,sites,pi_HA,pi_LA,dxy_HA_LA,Fst_HA_LA' tmp1.FST.popgenWindows > tmp2.FST.popgenWindows
awk -F "," '{print $1,$2,$3,$4,$9}' tmp2.FST.popgenWindows > tmp3.FST.popgenWindows


## 278299 regions
## 2,782 regions  --1%
## 13,915 regions --5%
## 27830 regions  --10%
sort -n -k5 -r tmp3.FST.popgenWindows | sed -n '1,2782p' > top1
sort -n -k5 -r tmp3.FST.popgenWindows | sed -n '1,13915p' > top5
sort -n -k5 -r tmp3.FST.popgenWindows | sed -n '1,27830p' > top10

sort -n -k1 -k2 -k3 top1 > F.top1
sort -n -k1 -k2 -k3 top5 > F.top5
sort -n -k1 -k2 -k3 top10 > F.top10

sed  's/ /\t/g' F.top1 > F.top1.bed
sed  's/ /\t/g' F.top5 > F.top5.bed
sed  's/ /\t/g' F.top10 > F.top10.bed

<SET_YOUR_PATH_TO_SOFTWARE>/software/bedtools2/bin/bedtools intersect -a <SET_YOUR_PATH>/F.top1.bed -b <SET_YOUR_PATH_TO_SOFTWARE>/snpEff/data/sheep/genes.gff -wa -wb > top1.annotation

<SET_YOUR_PATH_TO_SOFTWARE>/software/bedtools2/bin/bedtools intersect -a <SET_YOUR_PATH>/F.top5.bed -b <SET_YOUR_PATH_TO_SOFTWARE>/snpEff/data/sheep/genes.gff -wa -wb > top5.annotation

<SET_YOUR_PATH_TO_SOFTWARE>/bedtools2/bin/bedtools intersect -a <SET_YOUR_PATH>/F.top10.bed -b<SET_YOUR_PATH_TO_SOFTWARE>/snpEff/data/sheep/genes.gff -wa -wb > top10.annotation
