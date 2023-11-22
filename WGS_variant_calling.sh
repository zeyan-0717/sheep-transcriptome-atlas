## raw data to clean data
java -jar $Trimmomatic PE -phred33 \
$path/BFC2016518-01-MXS9_HJMTHALXX_L4_1.fq.gz $path/BFC2016518-01-MXS9_HJMTHALXX_L4_2.fq.gz \
$clean_fastq/C_MXS9_1.fq.gz $unpaired_fastq/L_MXS9_1.fq.gz $clean_fastq/C_MXS9_2.fq.gz $unpaired_fastq/L_MXS9_2.fq.gz \
ILLUMINACLIP:$adapter/TruSeq3-PE.fa:2:30:10:1:TRUE LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 -threads 4 MINLEN:50

## alignment
$bam/bwa mem -1 -t 10 -k 32 -M -R "@RG\tID:MXS9\tSM:MXS9\tLB:lib\tPL:illumina" \
$ref/Oar_rambouillet.fa $clean_fastq/C_MXS9_1.fq.gz $clean_fastq/C_MXS9_2.fq.gz | $samtools/samtools view -@ 10 -bS -F 4 -f 0x2 -t $ref/Oar_rambouillet.fa.fai -o $bamfile/MXS9.T1.bam
$gatk/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" SortSam --INPUT $bamfile/MXS9.T1.bam -OUTPUT $bamfile/MXS9.T2.bam --SORT_ORDER coordinate
$gatk/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" MarkDuplicates --INPUT $bamfile/MXS9.T2.bam -M $bamfile/MXS9.marked_dup_metrics.txt --OUTPUT $bamfile/MXS9.duplicate.bam


## bam to gvcf
$gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
HaplotypeCaller -R $ref/Oar_rambouillet.fa \
-I $bam/'MXS9.duplicate.bam \
--genotyping-mode DISCOVERY \
--input-prior 0.001 --input-prior 0.4995 -ERC GVCF -OVI true -L '$j' \
-O $output/chr'${j}'/'${i}.${j}'.g.vcf.gz


## Joint genotyping gvcf
$gatk-4.1.2.0/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" \
GenomicsDBImport \
--genomicsdb-workspace-path /public/home/casdao/kdylvfenghua/kdylvfenghua/Adaptation_Project/vcf/chr'${i}'/DB \
--batch-size 100 \
--intervals '${i}' \
--sample-name-map /public/home/casdao/kdylvfenghua/kdylvfenghua/Adaptation_Project/gvcf/gvcf_T/859_map.list \
--tmp-dir=/public/home/casdao/kdylvfenghua/kdylvfenghua/Adaptation_Project/vcf/chr'${i}'/tmp \
--reference /public/home/casdao/kdylvfenghua/kdylvfenghua/reference/Oar_rambouillet.fa \
--reader-threads 10'


## output vcf file
$gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" \
GenotypeGVCFs \
--reference /public/home/casdao/kdylvfenghua/kdylvfenghua/reference/Oar_rambouillet.fa \
-V gendb:/public/home/casdao/kdylvfenghua/kdylvfenghua/Adaptation_Project/vcf/chr'${i}'/DB \
-O /public/home/casdao/kdylvfenghua/kdylvfenghua/Adaptation_Project/vcf/chr'${i}'/raw_chr'${i}'.vcf.gz'

