nextflow.enable.dsl=2

/** 
 running Yang's version of reachtools to combine R1 and R2 fq files
 Input must be .fq.gz NOT .fastq.gz - we need a check for this
*/
process REACHTOOLS_COMBINE2 {

    //label 'mem40G' 

    publishDir "${outputdir}/nf_chosen_outputs",
		mode: "link", overwrite: true

    input:
    tuple val(sample_id), path(reads)
    val (outputdir)

    output:
        tuple val(sample_id), path ("*combined.fq.gz"), emit: comb
        path "*type_report.xls.gz", optional: true, emit: report

    script:
    """
    module load paired_tag   
    mkdir reachtools_${sample_id}_logs   
    reachtools combine2 $sample_id
    """
}

process REACHTOOLS_COMBINE3 {
    publishDir "${outputdir}/nf_chosen_outputs",
		mode: "link", overwrite: true

    input:
    tuple val(sample_id), path(reads)
    val (outputdir)

    output:
        tuple val(sample_id), path ("*combined3.fq.gz"), emit: comb
        path "*type_report.xls.gz", optional: true, emit: report

    script:
    """
    module load paired_tag   
    mkdir reachtools_${sample_id}_logs   
    reachtools combine3 $sample_id
    """
}

process REACHTOOLS_COMBINE_BULK {
    publishDir "${outputdir}/nf_chosen_outputs",
		mode: "link", overwrite: true

    input:
        tuple val(sample_id), path(reads)
        val (outputdir)

    output:
        tuple val(sample_id), path ("*combined.fq.gz"), emit: comb
        path "*type_report.xls.gz", optional: true, emit: report

    script:
    """
    module load paired_tag   
    mkdir reachtools_${sample_id}_logs   
    reachtools combinebulk $sample_id
    """
}

process REACHTOOLS_RMDUP2 {

    label 'bigMem'

    publishDir "${outputdir}/nf_chosen_outputs",
		mode: "link", overwrite: true

    input:
        tuple val(sample_id), path(mapped_reads)
        val (outputdir)

    output:
        tuple val(sample_id), path ("*rmdup.bam"), emit: reads

    script:
    """
    module load paired_tag 
    module load samtools    
    reachtools rmdup2 $mapped_reads
    """
}

// when we're using this later on we don't have the read name anymore
// don't think we actually need this any more
// process REACHTOOLS_RMDUP2_ONLY_BAM {

//     publishDir "${outputdir}/nf_chosen_outputs",
// 		mode: "link", overwrite: true

//     input:
//         path(mapped_reads)
//         val (outputdir)

//     output:
//         path ("*rmdup.bam"), emit: reads

//     script:
//     """
//     module load paired_tag 
//     module load samtools    
//     reachtools rmdup2 $mapped_reads
//     """
// }