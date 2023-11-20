##  set dir 
# set foler
workdir=<SET_YOUR_DIR>

for id in $(cat $workdir/script/joint_sample) # you need to prepare your "joint_sample" file
do
echo "**** ${id} joint-calling peaks start ****"
Genrich -t $workdir/genrich/${id}_1.sortn.last.bam,$workdir/genrich/${id}_2.sortn.last.bam -o $workdir/genrich/${id}_merge.raw.narrowPeak \
-j \
-p 0.01 \
-v && echo "**** ${id} joint-calling peaks finish ! ****"
# remove redudant peaks
grep -v "^PEKD" $workdir/genrich/${id}_merge.raw.narrowPeak | grep -v "^KZ" > $workdir/genrich/${id}_merge.narrowPeak
# rm $workdir/genrich/${id}_merge.narrowPeak
rm $workdir/genrich/${id}_1.sortn.last.bam $workdir/genrich/${id}_2.sortn.last.bam
awk '{print $1"\t"$2"\t"$3"\t"$4"\t+"}' $workdir/genrich/${id}_merge.narrowPeak > $workdir/bed/${id}_merge_peak.bed && echo "** ${id} merged bed file extract finish **" # extract *bed file 
done
