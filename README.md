# Nextflow Pipeline for Yang's scripts

A previous version of this pipeline can be found here: https://github.com/laurabiggins/nf_yang_wang 

 
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
```
                                                                    --- READS_QC_EXTRACT_BC2
                                                                    
                                                         


Yang has a set of scripts that he uses to process sequencing data with multiple barcodes. We are putting this into a NextFlow pipeline to try and streamline the usage of it.
