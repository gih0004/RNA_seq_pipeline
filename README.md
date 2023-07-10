# HTSEQ_counts
RNA seq pipeline designed for paired end reads using HTSEQ for producing files needed for visualizing transcriptomic data with DESEQ2. To use this script, put all fastq files within a single directory containing additionally your reference genome fna file and your reference gtf file.  
Fastq file names should be formated as follows: sample_1.fq.gz sample_2.fq.gz

This script requires three arguments and in specific order within the command line when submititng the script to HPC and shouuld look like: 
```ruby
$ htseq_counts <species name> <reference genome fna> <reference gtf> 
```
where species name can be either a common name for the species or scientific, should be one string in total.
where refernce genome fna is a fasta file that is used as the reference genome for sample alignment and indexing 
where reference gtf is a gene transfer format  file neccesary for generarting HTSEQ count file 

```ruby 
!/bin/bash
#Make sure to change all <pwd> with the current working directory where you have all fastq raw reads and your gft and reference genome 
# Format for Fastq raw read file names : <sample>_1.fq.gz <sample>_2.fq.gz

#SBATCH --job-name=HTSEQ_counts
#SBATCH --ntasks=10
#SBATCH --partition=bigmem2
#SBATCH --export=ALL
#SBATCH --array=1-49
#SBATCH --time=96:00:00
#SBATCH --error=/<pwd>/error.err
#SBATCH --output=/<pwd>/output.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<easley associated email> 

These SBATCH commands are condtions specific for Easley HPC, the only important and universal command from this block is the shebang line. if using easley, to submit the script:
$ sbatch HTSEQ_counts <species> <reference genome> <reference gtf> 
```


Step one is preprocessing of rna raw reads. This block takes the fastq rna raw reads and does three things:  
1) runs fastqc on the raw reads   
2) trims of adapters using fastp and an automatic adapter recognition option  
3) runs fastqc on the filtered AND adapterterless rna reads created by fastp 
Adapters can also be specified within the fastp command

### STEP 1: Run fastqc and Filter Raw reads
```ruby
module load fastqc
mkdir ./FASTQC
fastqc *.fq.gz  -o ./FASTQC/  


for fq1 in ./*_1.fq.gz 
do 
    base="${fq1%_1.fq.gz}"
   
    fastp -i $fq1 -I ${base}_2.fq.gz -o ${base}.filtered.1.fq.gz -O ${base}.filtered.2.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 -h ${base}_fastp.html -j ${base}_fastp.json


done
mkdir ./FASTQC_filtered

cp *filtered* ./FASTQC_filtered
cd ./FASTQC_filtered

fastqc *.filtered.1*  
fastqc *.filtered.2* 
cd ..


mkdir FASTQC.html
mkdir FASTQC_filtered.html

cp ./FASTQC/*.html ./FASTQC.html
cp ./FASTQC_filtered/*.html ./FASTQC_filtered.html

```
Results of Step 1: You should have three new directories:  
a. FASTQC = has original fastqc from rna raw reads  
b. FASTQC_filtered  = has fastqc results from adapterless and quality filtered reads   
c. FASTQC_filtered.html = contains all html files for viewing fastqc results in browser from the filtered and adapterless reads 




### STEP 2: Run HISAT2

Hisat firs creates the index based on a fna file and then aligns reads to the created index.   
Usage:  
hisat2-build [options]* <reference_in> <ht2_base>  
<reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves = {reference_fa}   
<ht2-base> The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base> = {species}  
```ruby 

module load hisat2
hisat2-build ${reference_fa} ${species}
echo "Finished index creation" >> Prrogress_ 
```
Note the Progress_ file, its purpose is to document where within the pipeline is the HPC at that time. To view you can do `less Progress_` on the comand line. This is the same for the execution of the rest of the script. 
After creating indices from genome, we can then run alignment of the samples against the indices recently created:  
```ruby
for fq1 in ./*_R1_001.fastq.gz #This should be changed into whatever you have last in sample names common between all samples 
do
echo "working with file $fq1" #only neccesary for when running through the terminal
base=$(basename $fq1 _R1_001.fastq.gz )
echo "base name is $base"

fq1=./${base}_R1_001.fastq.gz
fq2=./${base}_R2_001.fastq.gz

hisat2 -p 8 --dta -x ./sunflower/sunflower -1 $fq1 -2 $fq2 -S ${base}.sunflower.sam
done
echo "HISAT2 finished running!" # only use this line if running directly from terminal
#-dta is for downstream aplications such as Stringtie
#-p is for processors being used
#-S is for the output sam file
#-x is for the indices built and being used for alignment
```



### Step 3: converting SAM files to BAM files
```ruby
module load samtools
for i in /scratch/aubclsb0203/Project/HaIM/*.sam  #you can change the absolute path to relative path
do
samtools sort  ${i}  -o ${i}.sort.bam
done
```

To convert  gff to gtf file, which wont always be neccesary, use following command: 
```ruby
module load gffread/0.9.8
gffread genomic.gff -T -o genomic_sunflower.gtf
```


### STEP 4: Run featureCounts - StringTie
To load the module first specify source and then module to load:
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load stringtie/1.3.3
```ruby
for file in *.sunflower.sam.sort.bam
do
       tag=${file%.sunflower.sam.sort.bam}
stringtie -p 8 -G genomic_sunflower.gtf -o $tag.gtf -l sunflower  $tag.sunflower.sam.sort.bam
done
```
Using a wildcard for gtf files made from step 4, create a merged txt file of all gtf files generated from stringtie called mergelist.txt 
```ruby
ls *L002.gtf > mergelist.txt 
```

The command below takes the merged text file recently created and returns a proper merged gtf file called stringtie_merged_sunflower.gtf:
```ruby
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load stringtie/1.3.3
stringtie --merge -p 8 -G genomic_sunflower.gtf -o stringtie_merged_sunflower.gtf mergelist.txt
```

### STEP 5: Generating count table for ballgown:
To use ballgown, a count table from all the gtf files must be made. Stringtie will compare each sample agianst merged assembly to espimate transcript abundance
```ruby
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load stringtie/1.3.3

for file in *.bam;

    do tag=${file%.bam};

stringtie -e -B -p 8 -G stringtie_merged_sunflower.gtf -o /scratch/aubclsb0203/Project/HaIM/ballgown/$tag/$tag.gtf $tag.bam

done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```
Now you will have neccesary output files to use in ballgown analysis, such files are : 
1. e_data.ctab
2. e2t.ctab
3. i_data.ctab
4. i2t.ctab
5. t_data.ctab
You should have these 5 files within a directory called ballgown and within subdirectories based of the sample names 
