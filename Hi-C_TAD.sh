workdir=<SET_YOUR_DIR>
trimmomatic=<SET_YOUR_DIR/Trimmomatic-0.39>
HiC_Pro=<SET_YOUR_DIR/HiC-Pro-master/HiC-Pro-master/bin> 
cworld_dekker=<SET_YOUR_DIR/cworld-dekker/scripts>

# download one Hi-C data from SRA under accession number SRR19426890
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR194/090/SRR19426890/SRR19426890_1.fastq.gz . && mv SRR19426890_1.fastq.gz Tibetan_sheep_Hi-C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR194/090/SRR19426890/SRR19426890_2.fastq.gz . && mv SRR19426890_2.fastq.gz Tibetan_sheep_Hi-C_2.fastq.gz


### step 1: remove adapters
cd <SET_YOUR_PATH/HI_C>
 id=tibetan_hic
 echo "**** trimmomatic cut adapters ${id} start *****"
 java -jar $trimmomatic/trimmomatic-0.39.jar PE \
 -threads 20 \
 -phred33 \
 $workdir/data/${id}_1.fastq.gz $workdir/data/${id}_2.fastq.gz \
 $workdir/data/${id}_paired_1.fastq.gz $workdir/data/${id}_unpair_1.fastq.gz \
 $workdir/data/${id}_paired_2.fastq.gz $workdir/data/${id}_unpair_2.fastq.gz \
 ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && echo "**** trimmomatic cut adapters ${id} finished  ****"
 check quality with fastQC
 fastqc -t 20 $workdir/data/${id}_paired_1.fastq.gz $workdir/data/${id}_paired_2.fastq.gz -o fastQC && echo "**** ${id} fastqc finished !****"


### step 2: alignment
# active environment
 source ~/.bashrc
 conda activate hic # activate HiC-Pro

# build bowtie2 index
 bowtie2-build $workdir/ref/sheep_v1.0.fa $workdir/ref/sheep_v1.0 && echo "**** sheep bowtie2 index build ****"

# create enzyme site file (here use DpnII enzyme)
 $HiC_Pro/utils/digest_genome.py $workdir/ref/sheep.fa -r DpnII -o $workdir/ref/sheep_DpnII.bed && echo "**** sheep enzyme site create ****"

# change config file

# start running
 mv $workdir/data/${id}_paired_1.fastq.gz $workdir/data/${id}_R1.fastq.gz
 mv $workdir/data/${id}_paired_2.fastq.gz $workdir/data/${id}_R2.fastq.gz
 $HiC_Pro/HiC-Pro --input $workdir/data/ --output $workdir/hicpro_output --conf $workdir/config-hicpro.txt && echo "**** Hi-C running finished ****"


### step 3: HiCPeak (create Cooler file)
for i in $(cat ~/2workspace/HI_C/ref/sheep_chr.txt) 	 
do 
grep "^$i" tibetan_hic_40000_abs.bed|awk 'NR==1{printf "'$i'" "\t" $4 "\t"}END{print $4}' >> sheep.bed
done

a=($(cut -f2 sheep.bed))	
b=($(cut -f3 sheep.bed))
mkdir 40K 


for i in $(cat ~/2workspace/HI_C/ref/sheep_chr.txt)	
do
awk -v x=${a[i-1]} -v y=${b[i-1]} '{if(x<=$1&&$1<=y&&$2>=x&&$2<=y){print $0}}' tibetan_hic_40000.matrix |awk -v x=${a[i-1]} '{print $1-x "\t" $2-x "\t" $3}'  >./40K/"$i"_"$i".txt
done


#### HiC Peak ####
source ~/.bashrc 
conda activate hicpeaks

input_path=~/2workspace/HI_C/hicpro_output/hic_results/matrix/tibetan_hic/raw/40000/40K
output_path=~/2workspace/HI_C/tad

echo "**** cool file generate start ****"
~/anaconda3/envs/hicpeaks/bin/toCooler -O $output_path/sheep_40000.cool \
-d $input_path/datasets \
--chromsizes-file $input_path/sheep_fa.size \
--nproc 20 && echo "**** cool file generate finish ****"

conda deactivate


#### TADLib #### 
source ~/.bashrc 
conda activate TADLib

~/anaconda3/envs/TADLib/bin/hitad -O $output_path/TAD_40K.txt \
-d $output_path/meta_file \
-p 30 \
--logFile $output_path/hitad.log && echo "**** TAD calling finish ! ****"

date
