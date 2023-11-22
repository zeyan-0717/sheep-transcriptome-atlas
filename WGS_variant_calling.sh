## set variables
gatk=<SET_YOUR_PATH_TO_SOFTWARE>
bwa=<SET_YOUR_PATH_TO_SOFTWARE>
samtools=<SET_YOUR_PATH_TO_SOFTWARE>
Trimmomatic=<SET_YOUR_PATH_TO_SOFTWARE>

ref=<SET_YOUR_PATH_TO_REFERENCE_GENOME>
clean_fastq=<SET_YOUR_PATH>
unpaired_fastq=<SET_YOUR_PATH>
bamfile=<SET_YOUR_PATH>
output=<SET_YOUR_PATH>


## raw data to clean data
java -jar $Trimmomatic PE -phred33 \
$path/<YOUR_DATA>_1.fq.gz $path/<YOUR_DATA>_2.fq.gz \
$clean_fastq/<YOUR_DATA>_1.fq.gz $unpaired_fastq/<YOUR_DATA>_1.fq.gz $clean_fastq/<YOUR_DATA>_2.fq.gz $unpaired_fastq/<YOUR_DATA>_2.fq.gz \
ILLUMINACLIP:$adapter/TruSeq3-PE.fa:2:30:10:1:TRUE LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 -threads 4 MINLEN:50

## alignment
$bam/bwa mem -1 -t 10 -k 32 -M -R "@RG\tID:MXS9\tSM:MXS9\tLB:lib\tPL:illumina" \
$ref/Oar_rambouillet.fa $clean_fastq/<YOUR_DATA>_1.fq.gz $clean_fastq/<YOUR_DATA>_2.fq.gz | $samtools/samtools view -@ 10 -bS -F 4 -f 0x2 -t $ref/Oar_rambouillet.fa.fai -o $bamfile/<YOUR_DATA>.T1.bam
# sort coordinate
$gatk/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" SortSam --INPUT $bamfile/<YOUR_DATA>.T1.bam -OUTPUT $bamfile/<YOUR_DATA>.T2.bam --SORT_ORDER coordinate
# remove the duplicates
$gatk/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" MarkDuplicates --INPUT $bamfile/<YOUR_DATA>.T2.bam -M $bamfile/<YOUR_DATA>.marked_dup_metrics.txt --OUTPUT $bamfile/<YOUR_DATA>.duplicate.bam


## bam to gvcf
$gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
HaplotypeCaller -R $ref/Oar_rambouillet.fa \
-I $bam/<YOUR_DATA>.duplicate.bam \
--genotyping-mode DISCOVERY \
--input-prior 0.001 --input-prior 0.4995 -ERC GVCF -OVI true -L '$j' \
-O $output/chr'${j}'/'${i}.${j}'.g.vcf.gz


## Joint genotyping gvcf
$gatk-4.1.2.0/gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" \
GenomicsDBImport \
--genomicsdb-workspace-path <YOUR_DATA_TO_VCF>/chr'${i}'/DB \
--batch-size 100 \
--intervals '${i}' \
--sample-name-map <YOUR_DATA_TO_GVCF>/gvcf_T/<YOUR_SAMPLE>.list \
--tmp-dir=<YOUR_DATA_TO_VCF>/chr'${i}'/tmp \
--reference $ref/Oar_rambouillet.fa \
--reader-threads 10'


## output vcf file
$gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=10" \
GenotypeGVCFs \
--reference $ref/Oar_rambouillet.fa \
-V gendb:<YOUR_DATA_TO_VCF>/chr'${i}'/DB \
-O <YOUR_DATA_TO_VCF>/chr'${i}'/raw_chr'${i}'.vcf.gz'

