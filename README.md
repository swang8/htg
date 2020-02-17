# High throughput genotyping with sequencing
Shichen Wang, Feb 18 2020.


## setup
```bash
mkdir $SCRATCH/genotyping
cd $SCRATCH/genotyping

pwd
```

## Download tools and sample data
1. Tools
> wget --no-check-certificate 'https://doc.google.com/uc?export=download&id=1Y2Imq_wrMIq2ZydgBCdCW_ioC6t9WtAY' -O tools.tar
> wget http://bit.ly/38EXFCM -O tools.tar

2. Data
> wget http://bit.ly/2uRMfNl -O data.tar

## Extract files
```bash
tar xvf tools.tar
mv Tools pl_scripts ~/

tar xvf data.tar
ls
```
## Analysis steps:

### 1. Prepare the reference genome
```
module load picard/2.15.0-Java-1.8.0
module load SAMtools/1.1-intel-2015B
module load Bowtie2

ref="chr6_ref.fa"
f=$(perl -e '$f=shift; $f=~s/\.fa$//; print $f, "\n"' $ref)

bowtie2-build $ref  $ref
java -jar  /sw/eb/software/picard/2.15.0-Java-1.8.0/picard.jar  CreateSequenceDictionary R=$ref O=${ref}.dict
ln -s ${ref}.dict ${f}.dict
samtools faidx $ref
ln -s ${ref}.fai ${f}.fai

```

### 1. Quality control
```
module load Trimmomatic/0.38-Java-1.8.0
# sample_1
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 2 $SCRATCH/genotyping/fastsq/sample_1_R1.fq.gz $SCRATCH/genotyping/fastsq/sample_1_R2.fq.gz  $SCRATCH/genotyping/fastsq/QC/sample_1_F.fq.gz  $SCRATCH/genotyping/fastsq/QC/sample_1_FU.fq.gz  $SCRATCH/genotyping/fastsq/QC/sample_1_R.fq.gz $SCRATCH/genotyping/fastsq/QC/sample_1_RU.fq.gz  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# QC for all samples
for input1 in `dir fastq/*R1.fq.gz`; do
  sample_name=$(perl -e '$f=shift; print $1 if $f=~/fastq\/(sample_\d+)/' $input1)
  input2=${input1/R1/R2}
  
  output_R1=$(dirname $input1)/QC/${sample_name}_F.fq.gz
  output_R1_single=$(dirname $input1)/QC/${sample_name}_FS.fq.gz
  output_R2=$(dirname $input1)/QC/${sample_name}_R.fq.gz
  output_R2_single=$(dirname $input1)/QC/${sample_name}_RS.fq.gz

  cmd="java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 2 $input1 $input2 $output_R1 $output_R1_single $output_R2 $output_R2_single ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

  echo $cmd
  # execute the command
  `$cmd`
done

```

### 2. Alignment
```bash
module load Bowtie2
# check the parameters for bowtie2
bowtie2 --help  2>&1 | less
```
> bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]


```bash
# sample_1
bowtie2 -x chr6_ref.fa -1 fastq/QC/sample_1_F.fq.gz -2 fastq/QC/sample_1_R.fq.gz --rg-id sample_1 --very-sensitive-local -p 10 -S sample_1.sam
```

### 3. Processing steps
1. Sam to bam, samtools
> samtools view -Sb sample_1.sam >sample_1.bam

2. sort bam, samtools
> samtools sort sample_1.bam -o sample_1_sorted.bam

3. Add RG name if necessary, Picard
> java -jar $picard_dir/AddOrReplaceReadGroups.jar I=sample_1_sorted.bam O=sample_1_sorted_addRG.bam PL=illumina PU=barcode SM=sample_1 LB=SeqCap ID=sample_1 VALIDATION_STRINGENCY=SILENT

4. Mark duplciate, Picard
> java -Xmx15G -Djava.io.tmpdir=./  -jar $picard_dir/MarkDuplicates.jar I=$bam O=$out M=$merics_file REMOVE_DUPLICATES=true AS=true

5. local realignment, GATK
> java -Xmx15G -jar $GATK_jar  -I $_ -R $ref_fasta -T RealignerTargetCreator -o $interval
> java -Xmx15G -jar $GATK_jar -I $_ -R $ref_fasta -T IndelRealigner -targetIntervals $interval -o $realn_bam

6. Mapping Quality filtering, samtools
> samtools view -h -q $maq $input

### 4. Calling variations with GATK
```bash
# unifiedgenotyper
java -Xmx30G  -jar $GATK_jar -T UnifiedGenotyper -R $ref_fasta -I sample_1_sorted_addRg_realign_MQ20.bam --genotype_likelihoods_model BOTH  -o Variations/raw -rf BadMate -rf DuplicateRead -U ALLOW_N_CIGAR_READS

# haplotypecaller
java -Xmx30G  -jar $GATK_jar -T HaplotypeCaller -R $ref_fasta -I sample_1_sorted_addRg_realign_MQ20.bam --genotype_likelihoods_model BOTH  -o Variations/raw -rf BadMate -rf DuplicateRead -U ALLOW_N_CIGAR_READS

```

<hr>

## Autoamte the whole process
1. Generate bsub files for all the steps
We will need to have 



