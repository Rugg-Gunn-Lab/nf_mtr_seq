# Nextflow Pipeline for processing scMTR-seq (single-cell Multi-Targets and mRNA sequencing) data

## Description

This pipeline processes single cell sequencing data generated using the method described in [Combinatorial profiling of multiple histone modifications and transcriptome in single cells using scMTR-seq](https://doi.org/).
A detailed bench protocol for performing this technique can be found here: [protocols.io doi](https://doi.org/)

Profiling combinations of histone modifications identifies gene regulatory elements in different states and discovers features controlling transcriptional and epigenetic programmes.
However, efforts to map chromatin states in complex,
heterogeneous samples are hindered by the lack of methods that can profile multiple histone modifications together with transcriptomes inindividual cells.

scMTR-seq (single cell multi-targets and mRNA sequencing) is a multi-omics sequencing technology that can simultaneously profile multiple histone modifications or other chromatin-bound proteins and the transcriptome in the same single cells.
scMTR-seq has high sensitivity, high cell recovery and is highly scalable in terms of starting material requirement.
scMTR-seq has high signal specificity when compared with that of other methods which whilst they may have good on traget signal detection had higher levels of off-target signal, see our manuscript for additional details.

The protocol used to produce the libraries analysed by this pipeline follow these high level steps:

i) pre-assemble antibodies specific for each target histone modification or protein with indexed proteinA-Tn5-adapters
ii) perform in situ Tn5-mediated tagmentation with indexed complexes
iii) capture nuclear mRNA with a barcoded poly-T primer followed by in situ reverse transcription
iv) label single cells using split-pool combinatorial barcoding

scMTR-seq has been applied to profile up to 6 histone modifications and transcriptomes in single cells from a broad range of heterogeneous samples including pluripotent stem cell differentiation and mouse pre-implantation embryos.

Here is an outline of the steps performed by this pipeline:

- Cell barcodes (BC) and Unique Molecular Identifiers (UMIs) are extracted from reads with a [modified version](https://doi.org/10.5281/zenodo.13833575) of [reachtools](https://github.com/cxzhu/Paired-Tag)
	- Sequences are extracted from fixed positions in Read2:
		- 1-8 BC3, 39-46 BC2, 77-84 BC1
		- 115-119
			- RT-ID in RNA data (Reverse Transcript Index)
			- antibody-ID in DNA (Index identifying the antibody)
		- 120-127 UMI 
	  Combinations of the complete barcode sequences are of the form: `BC1:BC2:BC3(RT-ID|antibody-ID):UMI`
- Extracted barcode sequences were mapped to all possible barcode combinations with Bowtie.
- These identifiers extracted from Read2 are stored in the names of the sequences in the fastq files (see example below) !! example to be included
- `Trim-galore` is use to trim adapter sequences with default settings
- Trimmed reads are mapped to the specified genome using `Bowtie2` for DNA and `STAR` for RNA with default settings
- `samtools` is used to filter for uniquely mapping reads of MAPQ>50 for RNA and MAPQ>10 for DNA and to retain only uniquely mapping reads with the correct orientation.
- Deduplication was performed based on the mapping positions of UMI and cell ID for RNA data and cell ID only for DNA
- Deduplicated RNA were featurecounted with a custom `splitpoolquantitation` script and converted to a sparse matrix and deduplicated DNA were converted to fragment files.

### Example Barcode Information Extraction

Raw Fastq file from a DNA library

Read1

```
@AV240405:AV_B_YW6191_PE150_27022025:2428602470:1:10102:0226:0064 1:N:0:1 ACTTGA_GTCGGACT
TTGTAGACCAGGCTGGACTCAAACTCAGAGATCTACCTGTGTCTGCCTCCCTCAGTACTGAAATTAAAGGCGTGAGCCACCACGTCCAGCTAAACTCAG
+
GLLLLKLLMMMMMMMLMMNNNNNNNNNNNNNNNNNNNMNMNNLNNNNMMMMMMNNNMMMNNMMNNMNMMMMMMMMMMMMMMMMMMLMMMMMMMMMMMML
```
Read2

```
@AV240405:AV_B_YW6191_PE150_27022025:2428602470:1:10102:0226:0064 2:N:0:1 ACTTGA_GTCGGACT
CCAGTTCAGTGGCCGATGTTTCGCACGGCGTACGACTTCTTCACAGGATTCGAGGAGCGTGTGCGAACTCAGACCCCTCTATCATCCACGTGCTTGAGAGGCCAGAGCATTCGCCATACTCCCAATAGATGTGTATAAGAGACAGACCTGAGGGTGTGAGCAAAAATCTGGAAAGTCAGTGGTGGTAAAATCTGGTACA
+
FGGGHGGIMMMM@MMMKMMMMMMKGNJHNNNMMNININMLMMMMMKLNNFLLNLKMKMNMNNMMILMMN>JMLM>EHMJCKMCIM>MM=MMH>C>JLMKLKMCCLGLE1GMLM=K=GIJGLKK@JKHC5IJFIJ9CIIJFIDJIHGFII=JHJG6JJJGI9/JJHJ39ICIBJIAFI0JBJJ@3J6JJ;IH1@AJGIJI
```

Identifier positions:

```
                                                                                                            Antibody-
| BC3  |                              | BC2  |                             | BC1  |                              | ID|| UMI  |                    | Genomic Sequence                                  |
CCAGTTCAGTGGCCGATGTTTCGCACGGCGTACGACTTCTTCACAGGATTCGAGGAGCGTGTGCGAACTCAGACCCCTCTATCATCCACGTGCTTGAGAGGCCAGAGCATTCGCCATACTCCCAATAGATGTGTATAAGAGACAGACCTGAGGGTGTGAGCAAAAATCTGGAAAGTCAGTGGTGGTAAAATCTGGTACA
^      ^                              ^      ^                             ^      ^                              ^   ^       ^                    ^                                                   ^
         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         
         10        20        30        40        50        60        70        80        90        100       110       120       130       140       150       160       170       180       190       
```

Combine Read1 and Read2, extract barcode for mapping:

```
@AV240405:AV_B_YW6191_PE150_27022025:2428602470:1:10102:0226:0064:CTCCCAAT:TTGTAGACCAGGCTGGACTCAAACTCAGAGATCTACCTGTGTCTGCCTCCCTCAGTACTGAAATTAAAGGCGTGAGCCACCACGTCCAGCTAAACTCAG:GLLLLKLLMMMMMMMLMMNNNNNNNNNNNNNNNNNNNMNMNNLNNNNMMMMMMNNNMMMNNMMNNMNMMMMMMMMMMMMMMMMMMLMMMMMMMMMMMML:CCTGAGGGTGTGAGCAAAAATCTGGAAAGTCAGTGGTGGTAAAATCTGGTACA:FII=JHJG6JJJGI9/JJHJ39ICIBJIAFI0JBJJ@3J6JJ;IH1@AJGIJI
CCTCTATCTCTTCACACCAGTTCACCATA
+
FGGGHGGIMLMMMMMKHMJCKMCIK=GIJ
```

Identifier Positions:

```
| BC3  || BC2  || BC1  ||UMI|
CCTCTATCTCTTCACACCAGTTCACCATA
^      ^^      ^^      ^^   ^
         |         |         
         10        20        
```

After barcode mapping, now ready for genome mapping with both R1 and R2 (before trimming)

Read1

```
@AV240405.1.10102.0226.0064:46:32:17:03:CTCCCAAT
TTGTAGACCAGGCTGGACTCAAACTCAGAGATCTACCTGTGTCTGCCTCCCTCAGTACTGAAATTAAAGGCGTGAGCCACCACGTCCAGCTAAACTCAG
+
GLLLLKLLMMMMMMMLMMNNNNNNNNNNNNNNNNNNNMNMNNLNNNNMMMMMMNNNMMMNNMMNNMNMMMMMMMMMMMMMMMMMMLMMMMMMMMMMMML
```

Read2

```
@AV240405.1.10102.0226.0064:46:32:17:03:CTCCCAAT
CCTGAGGGTGTGAGCAAAAATCTGGAAAGTCAGTGGTGGTAAAATCTGGTACA
+
FII=JHJG6JJJGI9/JJHJ39ICIBJIAFI0JBJJ@3J6JJ;IH1@AJGIJI
```

Usage with test data:
```
./nf_mtr_seq -bg --type dna --genome GRCh38 --condition bulk      --outdir test_bulk_results      data/test.HIST_Bulk_D0*fq.gz
./nf_mtr_seq -bg --type dna --genome GRCh38 --condition paired    --outdir test_paired_results    data/test.HIST_D_S1_*.fq.gz
./nf_mtr_seq -bg --type dna --genome GRCh38 --condition stringent --outdir test_stringent_results data/test.HIST_D_S1_*.fq.gz
./nf_mtr_seq -bg --type dna --genome GRCh38 --condition relax     --outdir test_relax_results     data/test.HIST_HIST_D_S1_*.fq.gz

./nf_mtr_seq -bg --type rna --genome GRCh38 --condition stringent --outdir test_rna data/test.SLX23062_AGTTCC_GTAAGGAG_HIST_R_S1_S1_L001_R*.fq.gz
```

Dependencies (versions listed are those used during development of the pipeline):

- nextflow (v23.10.1)  
- samtools (v1.19.2)  
- bowtie (v1.3.1)    
- bowtie2 (v2.5.3)  
- star (2.7.11b)  
- fastqc (v0.12.1)  
- trim_galore (v0.6.10)  
- python (v3.12.2)  
- bedtools (v2.31.0)   
- perl (v5.32.1)  
- reachtools https://zenodo.org/doi/10.5281/zenodo.13833575  
 
Here is a graphical representation of the workflow:
```
--- REACHTOOLS_COMBINE
    |
    --- TRIM_GALORE (if condition == "stringent" or "paired" ) else straight to BOWTIE_BC
        |
        --- BOWTIE_BC - Map barcodes with bowtie
            |
            --- SAM_TO_FASTQ or SAM_TO_COV_FASTQ_PAIRED to convert mapped barcodes back to fastq
                |
                --- TRIM_GALORE
                    |
                        subworkflow
                    --- { RNA }
                        | 
                        --- STAR
                            | 
                            --- RNA_MAP_PROCESS (sub workflow) 
                            --- REACHTOOLS_RMDUP2
                                |
                                --- SAMTOOLS_INDEX
                                --- ADD_BAM_CB_TAG
                                    |
                                    --- BAM2COUNT_MATRIX/BAM2COUNT_MATRIX_SPARSE
                                --- BAM2BED
                                    |
                                    --- READS_QC_EXTRACT_BC
                                --- SUMMARIZE_MAPPED_READ_CELLS     
                    --- { DNA }
                        |
                        --- BOWTIE2
                            | 
                            --- SAMTOOLS SORT (+ filter if condition == "stringent")
                                |
                                    subworkflows
                                --- { DNA_PAIRED_BULK } (if condition == "paired" or "bulk" )
                                    |
                                    --- SUMMARIZE_MAPPED_READ_CELLS
                                    --- REACHTOOLS_RMDUP2
                                        |
                                        --- SAMTOOLS_SORT_FILT
                                            |
                                            --- BAM2BED
                                                |
                                                --- FILTER_READS
                                                    |
                                                    --- READS_QC_EXTRACT_BC
                                                        |
                                                        --- SORT_DEDUP

                                --- DNA_STRINGENT_RELAX (if condition == "stringent" or "relax")
                                    |
                                    --- SUMMARIZE_MAPPED_READ_CELLS
                                    --- REACHTOOLS_RMDUP2
                                        |
                                        --- SAMTOOLS INDEX
                                            |
                                            --- BAM2BED
                                                |
                                                --- READS_QC_EXTRACT_BC
                                                    |
                                                    --- COUNT_PILEUPS
                                                        |
                                                        --- REMOVE_PILEUPS
                                                            |
                                                            --- SAMTOOLS_INDEX2
                                                                |
                                                                --- BAM2BED2
                                                                    |
                                                                    --- READS_QC_EXTRACT_BC2
```
                                                                    

A previous version of this pipeline can be found here: https://github.com/laurabiggins/nf_yang_wang 

**Citation:**
```
Combinatorial profiling of multiple histone modifications and transcriptome in single cells using scMTR-seq.
Yang Wang, Jingyu Li, Andrew A. Malcolm, William Mansfield, Stephen J. Clark, Ricard Argelaguet, Laura Biggins, Richard Acton, Simon Andrews, Wolf Reik, Gavin Kelsey, Peter J. Rugg-Gunn
```

This project is released under the [GPL v3](./LICENSE) license, for licensing information of the individual tools used in this pipeline please refer to their own documentation.
