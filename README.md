# Sagerstrom_zebrafish_hindbrain

## 1. single cell RNAseq/ATACseq analyses for the manuscript:

## 2. CellRanger pipeline
Danio rerio GRCz11.99 fasta and gtf files downloaded from https://uswest.ensembl.org/Danio_rerio/Info/Index.
GFP sequence added to end of fasta file:
>chrGFP 
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCAGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCCACCGGATCTAGATAACTGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAGAAAACCTCCCACACCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCAAAAAAAAAAAAAAA

and end of gtf file:
chrGFP	unknown	gene	1	717	.	+	.	gene_id "GFP"; gene_version "2"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; 
chrGFP	unkonwn	transcript	1	717	.	+	.	gene_id "GFP"; gene_version "2"; transcript_id "GFP"; transcript_version "1"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; transcript_name "GFP"; transcript_source "manual"; transcript_biotype "protein_coding";
chrGFP	unknown	exon	1	717	.	+	.	gene_id "GFP"; gene_version "2"; transcript_id "GFP"; gene_name "GFP"; gene_source "manual"; gene_biotype "protein_coding"; transcript_name "GFP"; transcript_biotype "protein_coding";

To prevent conflict with chromosome naming in BSgenome.Drerio.UCSC.danRer11 used in Signac processing (ensembl used 1,2,3,etc... and BGgenome uses chr1, chr2, chr3, et...) added chr to the start of each line in ensembl gtf file and fasta label
