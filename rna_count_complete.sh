#!/bin/bash
#SBATCH -J RNA_quantify #job name
#SBATCH -p xhacnormala 
#SBATCH -N 1 # node number
#SBATCH -n 8
#SBATCH --cpus-per-task=4
#SBATCH -o quantify_report.out
date

## set path
# folder
workdir=<SET_YOUR_PATH>
index=<YOUR_PATH_TO_STAR_INDEX/star_index/>
# file
reference=<YOUR_PATH_TO_REFERENCE/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa>
gtf=<YOUR_PATH_TO_GTF/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf>

# step1:build index
STAR \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir $workdir/ref/star_index \
--genomeFastaFiles $reference \
--sjdbGTFfile $gtf \
&& echo "****index build****"

# step2: mapping to the genome
for id in $(cat $workdir/script/sample00) # you may need to change sample.txt
do
echo "**** STAR ${id} align start ****"
STAR \
--runThreadN 8 \
--runMode alignReads \
--readFilesCommand zcat \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFilterMismatchNmax 3 \
--outFilterMultimapNmax 10 \
--genomeDir $index \
--readFilesIn $workdir/clean_data/${id}_1.clean.fq.gz $workdir/clean_data/${id}_2.clean.fq.gz \ 
# check the suffix of the sequencing data
--outFileNamePrefix $workdir/align/${id}. \
&& echo "**** ${id} mapping finished ****"
# --quantMode TranscriptomeSAM GeneCounts

# step3: extract mapped reads
samtools view -h -f 2 -Sb $workdir/align/${id}.Aligned.sortedByCoord.out.bam > $workdir/align/${id}.bam
 # remove unimportant results
rm -rf $workdir/align/${id}._STARtmp 
rm $workdir/align/${id}.Log.progress.out $workdir/align/${id}.Log.out
# rm $workdir/align/${id}.Aligned.toTranscriptome.out.bam  # $workdir/align/${id}.Aligned.sortedByCoord.out.bam 

# step4: get read counts 
echo "**** quantify ${id} start ****"
featureCounts -T 7 -p -t exon -g gene_id -a $gtf -o $workdir/count/${id}.counts.txt $workdir/align/${id}.bam
# rm $workdir/count/${id}.counts.txt.summary

done


# step5: combine reports for QC ####
# results generate in "multiqc_data" folder under the corresponding current path
## summarize QC reports from fastqc
multiqc $workdir/fastqc/*.zip -o $workdir/fastqc/
## summarize QC reports from alignment software STAR
multiqc $workdir/align/*.final.out -o $workdir/align/ 
