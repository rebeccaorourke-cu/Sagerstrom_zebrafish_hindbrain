# Sagerstrom_zebrafish_hindbrain

## 1. single cell RNAseq/ATACseq analyses for the manuscript:
Kim YI, O'Rourke R, Sagerström CG. scMultiome analysis identifies embryonic hindbrain progenitors with mixed rhombomere identities. Elife. 2023 Nov 10;12. doi: 10.7554/eLife.87772. PubMed PMID: 37947350; PubMed Central PMCID: PMC10662952. https://pubmed.ncbi.nlm.nih.gov/37947350/

Data is available on GEO GSE223535 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223535

## 2. CellRanger pipeline
### 2.1 Make zebrafish GFP cellrannger reference
#### 2.1.1 Download reference files
Danio rerio GRCz11.99 fasta and gtf files downloaded from https://uswest.ensembl.org/Danio_rerio/Info/Index.

#### 2.1.2 Add GFP
GFP sequence added to end of fasta file:

>chrGFP 
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCAGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCCACCGGATCTAGATAACTGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAGAAAACCTCCCACACCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCAAAAAAAAAAAAAAA

and end of gtf file:

 chrGFP	unknown	gene	1	717	.	+	.	gene_id "GFP"; gene_version "2"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; 
 
 chrGFP	unkonwn	transcript	1	717	.	+	.	gene_id "GFP"; gene_version "2"; transcript_id "GFP"; transcript_version "1"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; transcript_name "GFP"; transcript_source "manual"; transcript_biotype "protein_coding";
 
 chrGFP	unknown	exon	1	717	.	+	.	gene_id "GFP"; gene_version "2"; transcript_id "GFP"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; transcript_name "GFP"; transcript_biotype "protein_coding";

#### 2.1.3 Correct chromosome naming
To prevent conflict with chromosome naming in BSgenome.Drerio.UCSC.danRer11 used in Signac processing (ensembl uses 1,2,3,etc... and BGgenome uses chr1, chr2, chr3, et...) added chr to the start of each line in ensembl gtf file and fasta ID line.

#### 2.1.4 Ran cellranger mkref

### 2.2 Cellranger arc 
Used cellranger_arc/1.0.1 count and zebrafishPlusGFP reference to generate gene expression matrix and ATAC fragment files.

## 3. Analyses in R with Signac
Used Signac/1.2.1 to analyze single cell data.  All analyses in workspace/notebooks



