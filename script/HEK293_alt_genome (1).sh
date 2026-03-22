# HEK293
# https://www.ebi.ac.uk/ena/browser/view/PRJNA399704
mkdir HEK293_genome
cd HEK293_genome

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR596/002/SRR5963472/SRR5963472_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR596/002/SRR5963472/SRR5963472_2.fastq.gz

# Map to the genome hg38
hisat2 --no-spliced-alignment -x hg38 -1 SRR5963472_1.fastq.gz -2 SRR5963472_2.fastq.gz -S SRR5963472.sam

# Remove sequencing reads
rm *.fastq.gz

# Conver to BAM and remove SAMs
samtools view -S -b SRR5963472.sam > SRR5963472.bam
rm *.sam

# Get the genome from UCSC for follow up analyses
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Index and create genome dictionary
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa

# Add read group information (required for mark duplicates)
samtools addreplacerg -r "ID:SRR1519317\tSM:1\tPL:ILLUMINA\tLB:lib1" -o SRR5963472_g.bam SRR5963472.bam

# Deduplication
# The Spark version permits the amalysis without BAM sorting. See GATK documentation.
gatk MarkDuplicatesSpark -I SRR5963472_g.bam -M dedup_metrics_1.txt  -O SRR5963472_dedup.bam 

# SNP calling. -ERC indicates GVCF input.
gatk HaplotypeCaller -ERC GVCF -R hg38.fa -I SRR5963472_dedup.bam -O SRR5963472_dedup_variants_p2.vcf.gz -ploidy 2
gatk HaplotypeCaller -ERC GVCF -R hg38.fa -I SRR5963472_dedup.bam -O SRR5963472_dedup_variants_p3.vcf.gz -ploidy 3

# Remove BAM files
rm *.bam*

# Genotyping
gatk GenotypeGVCFs -R hg38.fa -V SRR5963472_dedup_variants_p2.vcf.gz -O SRR5963472_genotype_p2.vcf.gz
gatk GenotypeGVCFs -R hg38.fa -V SRR5963472_dedup_variants_p3.vcf.gz -O SRR5963472_genotype_p3.vcf.gz

# Retrieve SNPs only
gatk SelectVariants -R hg38.fa -V SRR5963472_genotype_p2.vcf.gz --select-type-to-include SNP -O HEK293_SNP_raw_p2.vcf.gz
gatk SelectVariants -R hg38.fa -V SRR5963472_genotype_p3.vcf.gz --select-type-to-include SNP -O HEK293_SNP_raw_p3.vcf.gz

# Filter variants
gatk VariantFiltration -R hg38.fa -V HEK293_SNP_raw_p2.vcf.gz -O HEK293_SNP_tofilter_p2.vcf.gz --filter-name "AB_filter" --filter-expression "AB < 0.4"
gatk VariantFiltration -R hg38.fa -V HEK293_SNP_raw_p3.vcf.gz -O HEK293_SNP_tofilter_p3.vcf.gz --filter-name "AB_filter" --filter-expression "AB < 0.4"

# Remove filtered variants
zcat HEK293_SNP_tofilter_p2.vcf.gz | grep -v -e 'AB_filter' | bgzip > HEK293_SNPs_filtered_p2.vcf.gz
zcat HEK293_SNP_tofilter_p3.vcf.gz | grep -v -e 'AB_filter' | bgzip > HEK293_SNPs_filtered_p3.vcf.gz

# Index VCFs
gatk IndexFeatureFile -I HEK293_SNPs_filtered_p2.vcf.gz
gatk IndexFeatureFile -I HEK293_SNPs_filtered_p3.vcf.gz
# p2: 47950 SNPs, 1.5 per 100kb
# p3: 48510 SNPs, 1.5 per 100kb

# Build the altenative genome
gatk FastaAlternateReferenceMaker -R hg38.fa -V HEK293_SNPs_filtered_p2.vcf.gz -O HEK293p2_hg38.temp.fa
gatk FastaAlternateReferenceMaker -R hg38.fa -V HEK293_SNPs_filtered_p3.vcf.gz -O HEK293p3_hg38.temp.fa

# Revert to original chromosome names
sed -E 's/^>([0-9]+\s)([^:]+):[0-9\-]+/>\2/' HEK293p2_hg38.temp.fa > HEK293p2_hg38.fa
rm HEK293p2_hg38.fa.temp
sed -E 's/^>([0-9]+\s)([^:]+):[0-9\-]+/>\2/' HEK293p3_hg38.temp.fa > HEK293p3_hg38.fa
rm HEK293p2_hg38.temp.fa HEK293p3_hg38.temp.fa

# Remove other temps, compress genome
rm hg38.fa* hg38.dict
rm *temp*
rm *_dedup_* *_tofilter* *_SNP_raw*.vcf.gz* *_variants*.vcf.gz*

gzip HEK293p2_hg38.fa
gzip HEK293p3_hg38.fa


# # Indexing
# samtools faidx HEK293p2_hg38.fa
# gatk CreateSequenceDictionary -R HEK293p2_hg38.fa
# hisat2-build HEK293p2_hg38.fa HEK293p2_hg38
# samtools faidx HEK293p3_hg38.fa
# gatk CreateSequenceDictionary -R HEK293p3_hg38.fa
# hisat2-build HEK293p3_hg38.fa HEK293p3_hg38

