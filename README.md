# HTSEQ_counts
RNA seq pipeline designed for paired end reads using HTSEQ for producing files needed for visualizing transcriptomic data with DESEQ2. 
This script Requires three arguments and in specific order within the command line when submititng the script to HPC and shouuld look like: 
```ruby
$ htseq_counts <species name> <reference genome fna> <reference gtf> 
```
where species name can be either a common name for the species or scientific, should be one string in total.
where refernce genome fna is a fasta file that is used as the reference genome for sample alignment and indexing 
where reference gtf is a gene transfer format  file neccesary for generarting HTSEQ count file 

```ruby 
!/bin/bash
#Make sure to change all <pwd> with the current working directory where you have all fastq raw reads and your gft and reference genome 
# Format for Fastq raw read file names : <sample>_1.fq.fz <sample>_2.fq.gz

#SBATCH --job-name=pipeline
#SBATCH --ntasks=10
#SBATCH --partition=bigmem2
#SBATCH --export=ALL
#SBATCH --array=1-49
#SBATCH --time=48:00:00
#SBATCH --error=/<pwd>/error.err
#SBATCH --output=/<pwd/output.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<easley associated email> 

These SBATCH commands are condtions specific for Easley HPC, the only important and universal command from this block is the shebang line. if using easley, to submit the script:
$ sbatch HTSEQ_counts <species> <reference genome> <reference gtf> 
```

First, Change working directory into directory where all raw reads are stored in
`cd /scratch/aubclsb0203/Project/HaIM`

### STEP 1: Run fastqc
```ruby
module load fastqc
fastqc *.fastq.gz  -o /scratch/aubclsb0203/Project/HaIM #all fastq files within a directory will run fastqc on, it outputs to whatever path is given in -o option
```

Did not run trimmomatic to trim reads with poor quality because it was not neccesary due to quality of reads being good

For cleaning up GFF files downloaded from internet, and ensuring they only contain identifier in header row which is neccesary for creating indices, use: 

```ruby
awk -F "." '{print $1}' GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna > GCFedited.fna
```



### STEP 2: Run HISAT2
To create the genome indices: 

```ruby 
module load hisat2/2.0.5
hisat2-build --large-index  -f GCFedited.fna sunflower
```
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
