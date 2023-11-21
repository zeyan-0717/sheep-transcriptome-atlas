for i in {1..27}
do
bcftools view /public/home/casdao/kdylvfenghua/kdylvfenghua/climate_change/vcf/RD.chr${i}.snp.vcf.gz -c 1 -S /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/sample.list --threads 4 -Oz > /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/TIB.chr${i}.snp.vcf.gz
done 

for i in {1..27}
do
bcftools index -t /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/TIB.chr${i}.snp.vcf.gz
done 


for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/software/gatk-4.2.3.0/gatk SelectVariants -restrict-alleles-to BIALLELIC -V /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/TIB.chr${i}.snp.vcf.gz -O /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f1.TIB.chr${i}.snp.vcf.gz
done 


# filtering out MAF (Minor allele frequency) < 0.01 and Missing rate>0.05

for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/vcftools/bin/vcftools --gzvcf /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f1.TIB.chr${i}.snp.vcf.gz --maf 0.01 --max-missing 0.95 --recode --out /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz
done


for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/vcftools/bin/vcftools --vcf /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf --weir-fst-pop /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/HA_sheep.list --weir-fst-pop /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/LA_sheep.list  --out HA_vs_LA.chr${i} --fst-window-size 10000 --fst-window-step 10000
done

for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/htslib-1.12/bin/bgzip /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf -@ 6
done 

for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/htslib-1.12/bin/tabix  /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf.gz 
done

for i in {1..27}
do
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/Python3.8/bin/python3.8 /public/home/casdao/kdylvfenghua/kdylvfenghua/software/genomics_general-master/VCF_processing/parseVCFs.py -i /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf.gz --threads 6 | bgzip --threads 6 > /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.geno.gz
done

for i in {1..27}
do
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/Python3.8/bin/python3.8 /public/home/casdao/kdylvfenghua/kdylvfenghua/software/genomics_general-master/popgenWindows.py --windSize 10000 --stepSize 10000 --popsFile /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/popsFile  --outFile /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/FST.chr${i}
done

for i in {1..27}
do 
bcftools view -e 'F_MISSING > 0.05 & MAF <0.01' --threads 6 /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f1.TIB.chr${i}.snp.vcf.gz -Oz -o /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f3.TIB.chr${i}.snp.vcf.gz
done 

for i in {1..27}
do
/public/home/2020024/bin/bin/bcftools stats f1.TIB.chr${i}.snp.vcf.gz > f1.TIB.chr${i}.stats
done

for i in {1..27}
do
cat -n f1.TIB.chr${i}.stats | grep 'number of SNPs: ' >> 1_27_stats
done

for i in {1..27}
do 
bcftools index --threads 6 /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f3.TIB.chr${i}.snp.vcf.gz
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/htslib-1.12/bin/tabix  /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f3.TIB.chr${i}.snp.vcf.gz
done 


for i in {1..27}
do
source ~/.bashrc
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/Python3.8/bin/python3.8 /public/home/casdao/kdylvfenghua/kdylvfenghua/software/genomics_general-master/VCF_processing/parseVCFs.py -i /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f3.TIB.chr${i}.snp.vcf.gz --threads 6 --include ${i} | bgzip --threads 6 > /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.geno.gz
done

for i in {1..27}
do
source ~/.bashrc
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/Python3.8/bin/python3.8 /public/home/casdao/kdylvfenghua/kdylvfenghua/software/genomics_general-master/popgenWindows.py -g /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.geno.gz -f phased -T 5 --windSize 10000 --stepSize 10000 -p HA GNM327,GNM330,GNM341,CJL285,CJL289,CJL300,SAB311,SAB319,SAB322,ZJZ225,ZJZ226,ZJZ242,ZLX154,ZLX160,ZLX161,CHA02,CHA05 -p LA HU62,HU65,HU69,HU63,SXW30,SXW19,SXW29,SXW16,WDS123,WDS14,WDS109,WDS112,SSS44,SSS45,SSS46,SSS32,HDW47,HDW51,HDW52,HDW48  --popsFile /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/popsFile  --outFile /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/FST.chr${i}
done

for i in {1..27}
do
zcat f2.TIB.chr${i}.geno.gz | wc -l
done


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

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top1.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top1.annotation

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top5.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top5.annotation

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top10.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top10.annotation


-------------------------------------------------------------------------------------------
## 按照site-by-sites


for i in {1..27}
do 
/public/home/casdao/kdylvfenghua/kdylvfenghua/bin/vcftools/bin/vcftools --gzvcf /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/f2.TIB.chr${i}.snp.vcf.gz.recode.vcf.gz --weir-fst-pop /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/HA_sheep.list --weir-fst-pop /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/LA_sheep.list  --out site.HA_vs_LA.chr${i} 
done


cat site.HA_vs_LA.chr1.weir.fst  site.HA_vs_LA.chr2.weir.fst  site.HA_vs_LA.chr3.weir.fst  site.HA_vs_LA.chr4.weir.fst  site.HA_vs_LA.chr5.weir.fst  site.HA_vs_LA.chr6.weir.fst  site.HA_vs_LA.chr7.weir.fst  site.HA_vs_LA.chr8.weir.fst  site.HA_vs_LA.chr9.weir.fst  site.HA_vs_LA.chr10.weir.fst  site.HA_vs_LA.chr11.weir.fst  site.HA_vs_LA.chr12.weir.fst  site.HA_vs_LA.chr13.weir.fst  site.HA_vs_LA.chr14.weir.fst  site.HA_vs_LA.chr15.weir.fst  site.HA_vs_LA.chr16.weir.fst  site.HA_vs_LA.chr17.weir.fst  site.HA_vs_LA.chr18.weir.fst  site.HA_vs_LA.chr19.weir.fst  site.HA_vs_LA.chr20.weir.fst  site.HA_vs_LA.chr21.weir.fst  site.HA_vs_LA.chr22.weir.fst  site.HA_vs_LA.chr23.weir.fst  site.HA_vs_LA.chr24.weir.fst  site.HA_vs_LA.chr25.weir.fst  site.HA_vs_LA.chr26.weir.fst  site.HA_vs_LA.chr27.weir.fst > FST.popgen.sites



## 32786641 sites
## 327866   sites  --1%
## 1639332  sites  --5%
## 3278664  sites  --10%
sort -g -k3 -r FST.popgen.sites | sed -n '1,327866p' > top1
sort -g -k3 -r FST.popgen.sites | sed -n '1,1639332p' > top5
sort -g -k3 -r FST.popgen.sites | sed -n '1,3278664p' > top10

## 0.210384  --1%
## 0.114736  --5% 
## 0.0733178 --10%

sort -n -k1 -k2 -k3 top1  > sites.F.top1
sort -n -k1 -k2 -k3 top5  > sites.F.top5
sort -n -k1 -k2 -k3 top10 > sites.F.top10

sed  's/ /\t/g' F.top1 > F.top1.bed
sed  's/ /\t/g' F.top5 > F.top5.bed
sed  's/ /\t/g' F.top10 > F.top10.bed

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top1.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top1.annotation

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top5.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top5.annotation

/public/home/casdao/kdylvfenghua/kdylvfenghua/software/bedtools2/bin/bedtools intersect -a /public/home/casdao/kdylvfenghua/kdylvfenghua/YZ/F.top10.bed -b /public/home/casdao/kdylvfenghua/kdylvfenghua/software/snpEff/data/sheep/genes.gff -wa -wb > top10.annotation


-------------------------------------------------------------------------------------------
## 按照ZOOM specific regions


/public/home/2020024/conda_lyu/bin/vcftools --gzvcf f2.TIB.chr3.snp.vcf.gz.recode.vcf.gz  --weir-fst-pop HA_sheep.list --weir-fst-pop LA_sheep.list  --out chr3 #FST
/public/home/2020024/conda_lyu/bin/vcftools --gzvcf f2.TIB.chr3.snp.vcf.gz.recode.vcf.gz  --site-pi --keep HA_sheep.list  --out chr3.HA_sheep                   #HA_pi
/public/home/2020024/conda_lyu/bin/vcftools --gzvcf f2.TIB.chr3.snp.vcf.gz.recode.vcf.gz  --site-pi --keep LA_sheep.list  --out chr3.LA_sheep                   #LA_pi

/public/home/2020024/conda_lyu/bin/vcftools --gzvcf f2.TIB.chr3.snp.vcf.gz.recode.vcf.gz  --freq --keep HA_sheep.list  --out chr3.HA_sheep_fr                  #HA_freq
/public/home/2020024/conda_lyu/bin/vcftools --gzvcf f2.TIB.chr3.snp.vcf.gz.recode.vcf.gz  --freq --keep LA_sheep.list  --out chr3.LA_sheep_fr                  #LA_freq


cat chr3.weir.fst | sed -n '2462466,2473887p' > APOLD1.fst ## start:216823957; end:217823770 

cat chr3.HA_sheep.sites.pi | sed -n '2462515,2473936p' > APOLD1.HA.pi
cat chr3.LA_sheep.sites.pi | sed -n '2462546,2473969p' > APOLD1.LA.pi


cat chr3.LA_sheep_fr.frq | grep -n "216823957"  # 2462595
cat chr3.LA_sheep_fr.frq | grep -n "217823770"  # 2474018

cat chr3.HA_sheep_fr.frq | grep -n "216823957"  # 2462595
cat chr3.HA_sheep_fr.frq | grep -n "217823770"  # 2474018

cat chr3.LA_sheep_fr.frq | sed -n '2462595,2474018p' > APOLD1.LA.frq
cat chr3.HA_sheep_fr.frq | sed -n '2462595,2474018p' > APOLD1.HA.frq







