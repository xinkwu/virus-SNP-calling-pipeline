
# Modify according to your own situation
# software used
srapath=PATH_TO_sra
workspace=PATH_TO_workspace_bin
fastp=PATH_TO_fastp
fastqdump=PATH_TO_fastq-dump
samtools=PATH_TO_samtools
nanofilt=PATH_TO_nanofilt
minimap2=PATH_TO_minimap2
bowtie2=PATH_TO_bowtie2
pigz=PATH_TO_pigz
python3=PATH_TO_py3
# edit_sequence_using_indel.py edit_sequence_using_snps.py两个脚本放在workspace目录下
# genome & index used
humanbowtie2Reference=PATH_TO_human_bowtie2reference
humangenome=PATH_TO_humangenome_fasta
virusbwaReference=PATH_TO_virus_bwaReference
virusgenome=PATH_TO_virus_fasta

mkdir $workspace/fastq
mkdir $workspace/fastp
mkdir $workspace/maphuman
mkdir $workspace/rmhostreads
mkdir $workspace/mapvirus
mkdir $workspace/vcffiles
# added
mkdir $workspace/editedfasta

# parameter passed
sraid=$1
threadUse=$2

###### SNP calling pipeline ######

### Scenario 1: single-end Illumina
# Step1: raw reads decompress
$fastqdump --split-3 --outdir $workspace/fastq ${srapath}/${sraid}.sra
# Step2: reads QC
mkdir $workspace/fastp/$sraid
cd $workspace/fastp/$sraid && $fastp -w $threadUse -i $workspace/fastq/${sraid}.fastq -o ./${sraid}.clean.fastq.gz
# Step3: remove reads aligned to the host genome
$bowtie2 -q -p $threadUse --un-gz $workspace/rmhostreads/$sraid.clean.gz -x $humanbowtie2Reference -U $fastqpath/$sraid.clean.fastq.gz -S $workspace/maphuman/$sraid.human.sam
# Step4: alignment
$minimap2 -ax sr $virusgenome $workspace/rmhostreads/$sraid.clean.gz|$samtools view -hbS -F 4|$samtools sort > $workspace/mapvirus/${sraid}.bam
# # Step5: vcf calling
$samtools mpileup -A -f $virusgenome -Q 30 -u -v -t ADF,ADR $workspace/mapvirus/${sraid}.bam > $workspace/vcffiles/${sraid}.bam.vcf


### Scenario 2: paired-end Illumina
# Step1: raw reads decompress
$fastqdump --split-3 --outdir $workspace/fastq ${srapath}/${sraid}.sra
# Step2: reads QC
mkdir $workspace/fastp/$sraid
cd $workspace/fastp/$sraid && $fastp -w $threadUse -i $workspace/fastq/${sraid}_1.fastq -I $workspace/fastq/${sraid}_2.fastq -o ./${sraid}.clean.R1.fastq.gz -O ./${sraid}.clean.R2.fastq.gz
# Step3: remove reads aligned to the host genome
$bowtie2 -q -p $threadUse -x $humanbowtie2Reference -1 $workspace/fastp/$sraid/${sraid}.clean.R1.fastq.gz -2 $workspace/fastp/$sraid/${sraid}.clean.R1.fastq.gz -S $workspace/maphuman/$sraid.human.sam --un-conc $workspace/rmhostreads/$sraid.clean.fq
# Step4: alignment
$minimap2 -ax sr $virusgenome $workspace/rmhostreads/$sraid.clean.1.fq $workspace/rmhostreads/$sraid.clean.2.fq|$samtools view -hbS -F 4|$samtools sort > $workspace/mapvirus/${sraid}.bam
# Step5: vcf calling
$samtools mpileup -A -f $virusgenome -Q 30 -u -v -t ADF,ADR $workspace/mapvirus/${sraid}.bam > $workspace/vcffiles/${sraid}.bam.vcf


### Scenario 3: nanopore reads
# Step1: raw reads decompress
$fastqdump --split-3 --outdir $workspace/fastq ${srapath}/${sraid}.sra
# Step2: reads QC
mkdir $workspace/fastp/$sraid
cd $workspace/fastp/$sraid && cat $workspace/fastq/${sraid}.fastq |$nanofilt -q 10|gzip > ${sraid}.clean.fastq.gz
# Step3: remove reads aligned to the host genome
$minimap2 -ax map-ont $humangenome $workspace/fastp/$sraid/${sraid}.clean.fastq.gz > $workspace/maphuman/$sraid.human.sam
$samtools fastq -n -f 4 $workspace/maphuman/$sraid.human.sam |$pigz > $workspace/rmhostreads/$sraid.clean.gz
# Step4: alignment
$minimap2 -ax map-ont $virusgenome $workspace/rmhostreads/$sraid.clean.gz |$samtools view -hbS -F 4|$samtools sort > $workspace/mapvirus/${sraid}.bam
# Step5: vcf calling
$samtools mpileup -A -f $virusgenome -Q 30 -u -v -t ADF,ADR $workspace/mapvirus/${sraid}.bam > $workspace/vcffiles/${sraid}.bam.vcf





##### VCF file analysis #####

# 处理单突变类型的情况
grep -v -P "#|DP=0|INDEL" $workspace/vcffiles/${sraid}.bam.vcf|grep -P "\t[A-Z],<"|awk 'OFS="\t"{print $1,$2-1,$2,$4,$5,$3,$8}'|perl -wpe 's/,<\*\>//g'|perl -wpe 's/DP=[0-9]+;I16=([0-9]+),([0-9]+),([0-9]+),([0-9]+),\S+/$1\t$2\t$3\t$4/g'|awk 'OFS="\t"{if($7+$8+$9+$10>0){print $1,$2,$3,$4,$5,".",$7+$8+$9+$10,$9+$10,($9+$10)/($7+$8+$9+$10)}}'|awk '$9>0' > $workspace/vcffiles/${sraid}.var1.bed
# 两种突变情况
grep -v -P "#|DP=0|INDEL" $workspace/vcffiles/${sraid}.bam.vcf|grep -P "\t[A-Z],[A-Z],<"|awk 'OFS="\t"{print $1,$2-1,$2,$4,$5,$10}'|perl -wpe 's/\t([A-Z]),([A-Z]),<\S+/\t$1\t$2/g'|perl -wpe 's/\S+:([0-9]+),([0-9]+),([0-9]+),([0-9]+):([0-9]+),([0-9]+),([0-9]+),([0-9]+)$/$1\t$2\t$3\t$5\t$6\t$7/g'|awk 'OFS="\t"{print $1,$2,$3,$4,$5,".",$7+$8+$9+$10+$11+$12,$8+$11 "\n" $1,$2,$3,$4,$6,".",$7+$8+$9+$10+$11+$12,$9+$12}'|awk 'OFS="\t"{print $0,$8/$7}' > $workspace/vcffiles/${sraid}.var2.bed
# 三种突变的情况
grep -v -P "#|DP=0|INDEL" $workspace/vcffiles/${sraid}.bam.vcf|grep -P "\t[A-Z],[A-Z],[A-Z]"|awk 'OFS="\t"{print $1,$2-1,$2,$4,$5,$10}'|perl -wpe 's/\t([A-Z]),([A-Z]),([A-Z])\t/\t$1\t$2\t$3\t/g'|perl -wpe 's/\S+:([0-9]+),([0-9]+),([0-9]+),([0-9]+):([0-9]+),([0-9]+),([0-9]+),([0-9]+)$/$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8/g'|awk 'OFS="\t"{print $1,$2,$3,$4,$5,".",$8+$9+$10+$11+$12+$13+$14+$15,$9+$13 "\n" $1,$2,$3,$4,$6,".",$8+$9+$10+$11+$12+$13+$14+$15,$10+$14 "\n" $1,$2,$3,$4,$7,".",$8+$9+$10+$11+$12+$13+$14+$15,$11+$15}'|awk 'OFS="\t"{print $0,$8/$7}' > $workspace/vcffiles/${sraid}.var3.bed
cat $workspace/vcffiles/${sraid}.var1.bed $workspace/vcffiles/${sraid}.var2.bed $workspace/vcffiles/${sraid}.var3.bed|bedtools sort > $workspace/vcffiles/${sraid}.var.bed

# edit sequence
$samtools depth $workspace/mapvirus/${sraid}.bam -a|awk '$3=="0"'> $workspace/mapvirus/${sraid}.gap
$python3 $workspace/edit_sequence_using_snps.py $workspace/vcffiles/${sraid}.var.bed $workspace/editedfasta/${sraid}.edited.fasta $virusgenome 0.5 $workspace/mapvirus/${sraid}.gap # AF>0.5来进行筛选

# 考虑indel的情况
grep -P "INDEL" $workspace/vcffiles/${sraid}.bam.vcf|grep -v '#'|awk 'OFS="\t"{print $1,$2-1,$2,$4,$5,$3,$8}'|perl -wpe 's/DP=[0-9]+;I16=([0-9]+),([0-9]+),([0-9]+),([0-9]+),\S+/\t$1\t$2\t$3\t$4/g'|awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11}'|awk 'OFS="\t"{if($7+$8+$9+$10>0){print $1,$2,$3,$4,$5,".",$7+$8+$9+$10,$9+$10,($9+$10)/($7+$8+$9+$10)}}'|grep -v ',' > $workspace/vcffiles/${sraid}.indel.bed
less $workspace/editedfasta/${sraid}.edited.fasta|bioawk -c fastx '{print ">BetaCoV\/Wuhan-Hu-1\/2019\|EPI_ISL_402125""\n"$seq }' > $workspace/editedfasta/${sraid}.edited.clean.fasta
$python3 $workspace/edit_sequence_using_indel.py $workspace/vcffiles/${sraid}.indel.bed $workspace/editedfasta/${sraid}.final.fasta $workspace/editedfasta/${sraid}.edited.clean.fasta 0.5 # AF>0.5来进行筛选

