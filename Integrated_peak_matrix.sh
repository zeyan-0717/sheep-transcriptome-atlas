workdir=<SET_YOUR_DIR>

## step 1 : get peak read ##
# combine all *merge.narrowPeak file (bedtools merge)
cat *merge.narrowPeak > all_merged.narrowPeak
sort -n -k 1 -k 2 -k3 all_merged.narrowPeak | grep -v "X" > all_merged_sorted.narrowPeak # sort by position
grep "X" all_merged.narrowPeak | sort -n -k 2 -k 3 >> all_merged_sorted.narrowPeak
less -S all_merged_sorted.narrowPeak | grep -v "^PEKD" | grep -v "^KZ" > final_merged.narrowPeak # remove contig
bedtools merge -i final_merged.narrowPeak > all_merged.bed
rm all_merged_sorted.narrowPeak 
# for each tissue (use narrowPeak bed file)
cd $workdir/bed
tissue="heart, liver" # for example
for id in ${tissue}
do
cat *${id}_merge_peak.bed > ${id}_merge_peak.bed
sort -n -k 1 -k 2 -k3 ${id}_merge_peak.bed | grep -v "X" > ${id}_merge_sorted_peak.bed # sort by position
grep "X" ${id}_merge_peak.bed | sort -n -k 2 -k 3 >> ${id}_merge_sorted_peak.bed
less -S ${id}_merge_sorted_peak.bed | grep -v "^PEKD" |  grep -v "^KZ" > ${id}_final_merge_sorted_peak.bed # remove contig
bedtools merge -i ${id}_final_merge_sorted_peak.bed > ${id}_all_merged.bed
rm ${id}_merge_sorted_peak.bed 
done

# get peak read (multiBamSummary BED-file)
# The user provides a BED file that contains all regions that should be considered for the coverage analysis.
# A common use is to compare data coverages between two different samples for a set of peak regions
echo "**** get peak read start ****"
multiBamSummary BED-file -p 10 \
--BED $workdir/genrich/all_merge.bed \
--bamfiles $(for i in $(cat $workdir/bamfile);do echo "${i}";done) \
--labels $(for i in $(cat $workdir/bamfile_name);do echo "${i}";done) \
--outFileName $workdir/deeptools/readCounts_atac.npz \
--minMappingQuality 30 \
--outRawCounts $workdir/deeptools/readCounts_atac.tab && echo "**** get peak read finish ****"

# step 3: convert *.bam to *.bw file
for id in $(cat $workdir/sample_id)
do
echo "**** ${id} deeptools bamcoverage start ****"
bamCoverage -p 10 \
-b $workdir/final_align/${id}.last.bam \
-o $workdir/bw/${id}.bw \
-bs 10 \
--centerReads \
--normalizeUsing RPKM && echo "**** ${id} deeptools bamcoverage finish ****"
done
date

# get peak RPKM ()
echo "**** get peak RPKM start ****"
multiBigwigSummary BED-file -p 10 \
--bwfiles $(for i in $(cat $workdir/bwfile);do echo "${i}";done) \
--BED $workdir/genrich/all_merge.bed \
--labels $(for i in $(cat $workdir/bamfile_name);do echo "${i}";done) \
--outFileName $workdir/deeptools/readRPKM_atac.npz \
--outRawCounts $workdir/deeptools/readRPKM_atac.tab && echo "**** get peak RPKM finish ****"

date
