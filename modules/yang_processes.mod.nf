nextflow.enable.dsl=2

//script_path="/bi/group/bioinf/Laura_B/nf_yang_wang/scripts/"
//splitpool_path = "/bi/group/bioinf/Laura_B/nf_yang_wang/splitpoolquantitation/"

script_path = workflow.projectDir.toString() + "/scripts/"

def get_BC_genome(data_type) {

    scriptDir = workflow.projectDir

    def filePrefix = scriptDir.toString() + "/cell_id/cell_id." + data_type
    def filePrefixTest = filePrefix + ".1.ebwt"
    def testFile = new File(filePrefixTest)
    if (!testFile.exists()) {
        println("\nFile >>$filePrefixTest<< does not exist. Something's gone wrong with finding the barcode genome...\n")       
    } else {
        return filePrefix
    } 
}

// get the length of R2 in the fastq file
process GET_R2_LENGTH {

     input:
        tuple val(name), path(reads)

	output:
		env (R2_len),    emit: read2_length 

    script:
        if (reads instanceof List) {
            R2_file = reads[1]
        } else {
            println ("\n Read 2 fastq file not found. '\n\n")
            exit 1 
        }

        """
        R2_len="\$(zcat ${R2_file} | head -n 2 | awk 'NR%4==2 {print length}')"
        """
}

process SAM_TO_COV_FASTQ_PAIRED {

    input:
        tuple val(name), path(mapped_reads)
        val(R2_length)
        val (outputdir)

	output:
	    tuple val(name), path ("*fq.gz"),  emit: reads

   // publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """
        perl ${script_path}/sam2covFastq.paired.pl ${mapped_reads} ${R2_length}  
        """
}

process SAM_TO_FASTQ {

    //label 'mem40G'

    input:
        tuple val(name), path(mapped_reads)
        val (outputdir)

	output:
	    tuple val(name), path ("*fq.gz"),  emit: reads
		path "*log",                       emit: stats 

   // publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """
        module load paired_tag || echo "no module found"
        reachtools convert2 ${mapped_reads} 1>reads_processed.log
        """
}

//|| echo "no module found"

process ADD_BAM_CB_TAG {

    label 'bigMem'

    input:
        tuple val(name), path(mapped_reads)
        path(bai) // an index file is required even if it's not used directly
        val (outputdir)

    output:
        tuple val(name), path ("*10X.bam"), emit: reads

    publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:

        """ 
        module load python || echo "no module found"
        python3 ${script_path}/add.bam.CB.tag.py -i ${mapped_reads} -p 30  
        """
} 


// output is one file if --output matrix
// we might just want two diferent processes
process BAM2COUNT_MATRIX {

    label 'bigMem'

    input:
        tuple val(name), path(mapped_reads)
        val (gtf_file)
        val (outputdir)

    output:
        tuple val(name), path ("${name}_splitpool_output"), emit: splitpool

    publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:

        """ 
        module load python || echo "no module found"
        python3 ${script_path}/quantitate_splitpool_scrna.py --output matrix ${gtf_file} ${mapped_reads} ${name}_splitpool_output
        
        """
} 

// output is 3 files in a folder
process BAM2COUNT_MATRIX_SPARSE {

    label 'bigMem'

    input:
        tuple val(name), path(mapped_reads)
        val (gtf_file)
        val (outputdir)

    output:
        path("${name}_splitpool_output/matrix.mtx.gz"), emit: matrix
        path("${name}_splitpool_output/barcodes.tsv.gz"), emit: barcodes
        path("${name}_splitpool_output/features.tsv.gz"), emit: features

    publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:

        """ 
        module load python || echo "no module found"
        python3 ${script_path}/quantitate_splitpool_scrna.py --output sparse ${gtf_file} ${mapped_reads} ${name}_splitpool_output
        
        """
} 

process BAM2BED {

    //label 'bigMem' // 20GB

    input:
        tuple val(name), path(mapped_reads)
        path(bai) // better if the bam file has been indexed
        val (outputdir)
        val (bam2bed_args)

    output:
        tuple val(name), path ("*.bed.gz"), emit: reads

    //publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """ 
        module load bedtools || echo "no module found"
        bedtools bamtobed -i ${mapped_reads} ${bam2bed_args} | gzip > ${mapped_reads}.bed.gz
        """
}

process BAM2BED_NOINDEX {

    //label 'bigMem' // 20GB

    input:
        tuple val(name), path(mapped_reads)
        val (outputdir)
        val (bam2bed_args)

    output:
        tuple val(name), path ("*.bed"), emit: reads

   // publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """ 
        module load bedtools || echo "no module found"
        bedtools bamtobed -i ${mapped_reads} ${bam2bed_args} > ${mapped_reads}.bed
        """
}


process READS_QC_EXTRACT_BC {

    //label 'bigMem' // 20GB

    input:
        tuple val(name), path(bed_file)
        val (outputdir)

    output:
        tuple val(name), path ("*.tsv.gz"), emit: reads

    publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:

        """     
        perl ${script_path}/reads.qc.extract.bc.pl ${bed_file}
        """
}

process SUMMARIZE_MAPPED_READ_CELLS {

    //label 'mem40G'
    label 'hugeMem' // 80GB

    input:
        tuple val(name), path(fq_file)
        tuple val(name), path(bam)
        val (outputdir)

    output:
        tuple val(name), path ("*summary.xls"), emit: qc

    publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:
        // we only want R1 
        if (fq_file instanceof List) {
            fq = fq_file[0]
        } else {
            fq = fq_file
        }

        """   
        module load samtools || echo "no module found"
        perl ${script_path}/summarize_mapped_read_cells.pl ${fq} ${bam}
        """   
}

process COUNT_PILEUPS {

    input:
        tuple val(name), path(mapped_reads)
        val (outputdir)

    output:
        path ("*txt"), emit: tally

   // publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """ 
        module load python || echo "no module found"
        python3 ${script_path}/count_pileups.py ${mapped_reads} ${name}_filtered_pos_strand_tally.txt       
        """
}

process REMOVE_PILEUPS {

    input:
        tuple val(name), path(mapped_reads)
        val (tally) // this is the output from COUNT_PILEUPS
        val (cutoff)
        path (samtools_index) // we need this accessible
        val (outputdir)

    output:
        tuple val(name), path ("*bam"), emit: rmpile
        

   // publishDir "${outputdir}/nf_chosen_outputs",
	//	mode: "link", overwrite: true

    script:

        """ 
        module load python || echo "no module found"
        python3 ${script_path}/remove_pileups.py ${mapped_reads} ${tally} ${name}_sorted_rmdup_rmpiles.bam ${cutoff}
        """
}
