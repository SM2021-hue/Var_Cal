# Variant Calling
#bcftools##

##GATK##
#Starting from raw reads (.fastq) to calling variants (.vcf) and their downstream analysis.
## Tools Required:
#* BWA 
#* Picard 
#* GATK-4 
#* SnpEFF 
#* Samtools 

## Tool installation
#Install the required tools in a conda environment using “variant_call.yml” file.

#But before that, you need to have conda installed in the system. 
#If it's already not installed. To install conda, from the command line:
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

#if your system has miniconda3 do this

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
!bash miniconda.sh -b -p $HOME/miniconda3
##################open new terminal##########################
wget https://raw.githubusercontent.com/sk-sahu/notebooks/master/variant_call.yml
#Once the conda is installed. Now let's install the tools. This will take time because all the tools will be installed at once.

conda env create -f variant_call.yml 
#Now, If you type `conda env list` you should be able to see an environment called - variant_call. Change to that environment by entering
conda activate variant_call

#Download reference genome
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#download query sequences
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925749/SRR925749_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925749/SRR925749_2.fastq.gz
#Start analysis Create reference genome index
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

picard CreateSequenceDictionary R= Homo_sapiens.GRCh38.dna.primary_assembly.fa O= Homo_sapiens.GRCh38.dna.primary_assembly.dict 
#step-1 Alignment-map to reference
bwa mem -M -t 16 Homo_sapiens.GRCh38.dna.primary_assembly.fa SRR925749_1.fastq SRR925749_2.fastq > SRR925749.sam
#step-2 Generate sorted bam file
picard SortSam INPUT=SRR925749.sam OUTPUT=SRR925749_sorted_reads.bam SORT_ORDER=coordinate
#step-3 Check alignment summary
picard CollectAlignmentSummaryMetrics R=Homo_sapiens.GRCh38.dna.primary_assembly.fa I=SRR925749_sorted_reads.bam O=alignment_metrics.txt

picard CollectInsertSizeMetrics INPUT=SRR925749_sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf

samtools depth -a SRR925749_sorted_reads.bam > depth_out.txt
#step-4 Mark duplicates
picard MarkDuplicates INPUT=SRR925749_sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=insert_metrics.txt
#step-5 Add read groups
picard AddOrReplaceReadGroups I=dedup_reads.bam O=final.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
#step-6 build bam index
picard BuildBamIndex INPUT=final.bam
#step-7 call variants
gatk HaplotypeCaller -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -I final.bam -O raw_variants.vcf
#step-8 Extract SNPs and INDELs
#For SNPs
gatk SelectVariants -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -V raw_variants.vcf -select-type SNP -O raw_snps.vcf
#For INDELs
gatk SelectVariants -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -V raw_variants.vcf -select-type INDEL -O raw_indels.vcf

############################################################################

##freebayes##
#Here we will demonstarte variant calling using freebayes
#he steps 1 to 6 will be the same as of mentioned above
#Environment
conda install freebayes==1.3.1 -c bioconda
#Step7: Call Variants
freebayes -f ref_genome.fa final.bam > raw_variants.vcf
freebayes -f Homo_sapiens.GRCh38.dna.primary_assembly.fa final.bam > raw_freebayes_variants.vcf
#Filter VCFs
vcftools --vcf raw_freebayes_variants.vcf --minQ 20 --recode --recode-INFO-all --out fb_filtred_variants.vcf
##bcf-tools
#This#Here we will demonstarte variant calling using bcftools
#he steps 1 to 6 will be the same as of mentioned above
#Environment 
conda install bcftools==1.9 -c bioconda
samtools mpileup -ugf $WORKDIR/refgenome/Homo_sapiens.GRCh38.dna.chromosome.X.fa  $WORKDIR/bam/read.filter.rmdup.bam  |bcftools call -vmO z -o $WORKDIR/vcf/read.bcftools.vcf.gz
#Annotation with known variants
#Download required refernece dbsnp database in vcf format (for human from here: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/)
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz


#Compress and index your raw_variants.vcf file to be used by next step.

bgzip fb_filtred_variants.vcf.recode.vcf
tabix fb_filtred_variants.vcf.recode.vcf.gz


#also index downloaded reference.vcf
tabix 00-All.vcf.gz

#Annotate using bcftools
bcftools annotate -c ID \
        -a 00-All.vcf.gz fb_filtred_variants.vcf.recode.vcf.gz \
        > fb_filtred_variants_annot.vcf


