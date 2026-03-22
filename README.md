# Identification of Novel miRNA Target Sites in HEK293 Genome

## Project Overview
This project investigate the genomic differences between the **HEK293 cell line** and the **Human Reference Genome (GRCh38/hg38)**. By leveraging high-throughput sequencing data, I identified Single Nucleotide Polymorphisms (SNPs) that specifically alter **microRNA (miRNA) target sites** within the **3' Untranslated Regions (3'UTRs)**.

## Key Findings
- **Variant Identification:** Identified ~48,000 SNPs in HEK293 using GATK4, with over 18,500 located in regulatory 3'UTR regions.
- **Novel Regulatory Sites:** Discovered approximately 490 miRNAs with unique target sites not present in the hg38 reference.
- **Ploidy Impact:** Analyzed the differences between diploid and triploid states of HEK293, highlighting unique miRNA interactions (e.g., hsa-miR-6844).

## Bioinformatics Pipeline
The analysis was performed using a professional-grade computational workflow:
1. **Alignment:** HISAT2 for mapping reads to GRCh38.
2. **Variant Calling:** GATK4 (HaplotypeCaller) for high-confidence SNP detection on an HPC cluster.
3. **Target Prediction:** SeedVicious for scanning canonical miRNA seed matches.
4. **Data Integration:** R (tidyverse, biomaRt) for 3'UTR retrieval and SNP filtering.
5. **Functional Analysis:** Gene Ontology (GO) and KEGG pathway enrichment analysis.

## Repository Structure
- `/scripts`: BASH and R scripts for the GATK pipeline and downstream analysis.

