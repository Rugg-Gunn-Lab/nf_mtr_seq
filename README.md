# Nextflow Pipeline for Yang's scripts

Yang has a set of scripts that he uses to process sequencing data with multiple barcodes. We are putting this into a NextFlow pipeline to try and streamline the usage of it.

## Yang's run_pipeline.sh script

There are $argv arguments to determine which section of the script is run.

## 1. prepare
This creates a folder structure. I'll leave this for now as nextflow creates files itself. May need to come back to this if we need a certain folder structure.

## 2. merge
Joins all *R1.fastq.gz files to one file, same with R2

I don't think the merging needs to be in the pipeline, at least not initially - let's start the nextflow pipeline with fastq files that are ready to be processed.
We can create a nice standalone script for merging files if we need to.
For now, just run:   
`cat *R1.fastq.gz > xx_R1_merged.fq.gz`  
We need `.fq.gz` instead of `.fastq.gz` for the next step. 

---
---

## Starting from here for the nextflow pipeline

---

This is being implemented in the file yang1.nf

## 3. PRE DNA/RNA 
This runs the gz.sh script from the scripts folder.  

### Combining and extracting barcode reads - reachtools

It uses `reachtools combine`  
reachtools has now been installed on the cluster so it can be run using
```
module load paired_tag
reachtools combine xxxxx   
```

```
reachtools --help
    combine/combine2/combine3       Combine Read1 with Read2, extract BC and UMI to a merged fastq file. Combine2 is used with 7-bp length barcodes, while combine3 is used with 8-bp length barcodes.
```

Different parameters/commands are used depending on the stringency required	
and whether the samples are bulk or not  
`--norc = no mapping to reverse complement`  

### trim_galore

This is run if `$conditions== "stringent"` or `paired`, not for `relax` or `bulk`

### bowtie

The output from reachtools and trim_galore (if run) is mapped with bowtie.

The output from this step is a SAM file of mapped barcodes.  
The only difference between DNA/RNA is the bowtie reference.

### convert sam to fastq

To convert the mapped barcode reads back to fastq:

`reachtools convert2 ${s}_BC.sam`

There are also the individual perl scripts that contain this code.


## 4. RUN_trim

Runs trim_galore

## 5. RUN RNA
Maps using STAR
then runs various scripts - there's quite a lot in there  
- star  
- samtools to filter out low quality mappings
- bam2.10xbam.sh which runs add.bam.CB.tag.py which "Add barcodes and umi from Paired-Tag readname to field CB and UMI"  
- bam2countmatrix.sh runs quantitate_splitpool_scrna.py - https://github.com/s-andrews/splitpoolquantitation
- bedtools bamtobed
- perl scripts to sort and summarise

## 6. RUN_bam2mtx
This option just runs the `bam2countmatrix.sh` script which is also within some other steps.

## 7. RUN DNA

If paired or bulk, runs `summarize_mapped_read_cells.pl` - I'm assuming to produce an extra log/info file, not sure what the script does exactly.

- maps using bowtie2
- samtools to filter out low quality mappings and sort

if paired or bulk - couldn't use rmdup2
    - bedtools bamtobed
    - awk command to keep the read pairs that are on the same chromosome
    - cut, sort to extract the fragment related columns and cellids, umi
    - reads.qc.extract.bc.pl
    - sort, uniq to get unique chr, start, end, cellid, umi

if not paired or bulk - uses `reachtools rmdup2`
count_pileups
reads.qc.extract.bc.pl

## 8. RUN SUMMARY

Runs `summarize_mapped_read_cells.pl`

---

This is where a lot of the scripts come from:  
https://github.com/cxzhu/Paired-Tag   
There doesn't seem to be a licence in the repo


---

## Current status

I can't find any usage of `reactools combinebulk` other than in Yang's gz.sh script.



## Test nextflow files

1. run `salmon index`
    This requires a transcriptome.fa file

2. Run `salmon quant` 
