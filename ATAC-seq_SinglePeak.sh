## set dir
workdir=<SET_YOUR_PATH>
ref=<YOUR_PATH_TO_BWA/bwa_index/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa>
gtf=<YOUR_PATH_TO_GTF/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf>

# create folders
mkdir {clean_data,fastqc,align,rmdup,final_align,stat,peaks,deeptools,bed}

# shell 1: fq_RawBam.sh, process threads = 12
## fastqc ## -- step 1
for id in $(cat $workdir/script/sample00) # suffix of your sequencing data
do
fastqc -t 12 $workdir/clean_data/${id}_paired_1.fastq.gz $workdir/clean_data/${id}_paired_2.fastq.gz -o $workdir/fastqc && echo "**** ${id} QC check finished! ****"
# multiqc $workdir/fastqc/*fastqc.zip -z -o $workdir/fastqc && echo "**** QC check finished! ****"

## align ## -- step 2 
bwa mem -t 12 $ref \
$workdir/align/${id}_paired_1.fastq.gz \
$workdir/align/${id}_paired_2.fastq.gz | samtools sort -O bam -@ 12 -o - > $workdir/align/${id}.raw.bam && echo "**** ${id} alignment finish ****"
# samtools 
samtools index $workdir/align/${id}.raw.bam
samtools flagstat $workdir/align/${id}.raw.bam -@ 12  > $workdir/stat/${id}.raw.stat
bedtools bamtobed -i $workdir/align/${id}.raw.bam > $workdir/bed/${id}.raw.bed && echo "**** ${id}.last.bed finish ****"

#### remove PCR duplicate ####
# rmdup bam
echo "**** remove ${id} duplicates started at $(date) *****"  
gatk MarkDuplicates \
-I $workdir/align/${id}.bam \
-O $workdir/rmdup/${id}.rmdup.bam \
-M $workdir/rmdup/${id}.rmdup.metrics.txt \
-ASO coordinate \
--REMOVE_DUPLICATES true && echo "****${id} rmdup finish at $(date)****" 
# rmdup bam -- statistic
samtools flagstat $workdir/rmdup/${id}.rmdup.bam -@ 8 > $workdir/stat/${id}.rmdup.stat

# filtration & sort -- final.bam (based on location)
samtools view -h -f 2 -q 30 $workdir/rmdup/${id}.rmdup.bam |grep -v "MT" |samtools sort -O bam -@ 8 -o - > $workdir/final_align/${id}.last.bam && echo "**** ${id}.last.bam finish ****"
# final.bam -- index & statistic
samtools index $workdir/final_align/${id}.last.bam && echo "**** ${id}.last.bam index finish ****"
samtools flagstat $workdir/final_align/${id}.last.bam -@ 8 > $workdir/stat/${id}.last.stat && echo "**** ${id}.last.stat finish ****"
# bedtools -- last.bed
bedtools bamtobed -i $workdir/final_align/${id}.last.bam > $workdir/bed/${id}.last.bed && echo "**** ${id}.last.bed finish ****"
# sorted by name
samtools sort -n -@ 8 -O bam  -o $workdir/genrich/${id}.sortn.last.bam -T $workdir/genrich/${id}_samtools.tmp $workdir/final_align/${id}.last.bam && echo "**** ${id} sort by name finish****"
# single peak calling
Genrich -t $workdir/genrich/${id}.sortn.last.bam -o $workdir/genrich/${id}.raw.narrowPeak -j -p 0.01 -v  && echo "**** ${id} peak calling finish ****"
# remove redudant peaks
grep -v "^PEKD" $workdir/genrich/${id}.raw.narrowPeak | grep -v "^KZ" > $workdir/genrich/${id}.narrowPeak
# peak bed file
awk '{OFS="\t"} {print $4,$1,$2,$3,"+"}' $workdir/genrich/${id}.narrowPeak > $workdir/bed/${id}.peak.bed
# rm $workdir/genrich/${id}.raw.narrowPeak

done

done
