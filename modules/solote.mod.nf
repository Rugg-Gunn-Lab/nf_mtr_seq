nextflow.enable.dsl=2

params.no_output = false

process SOLOTE {

    tag "$name"

    label 'hugeMem'
    label 'multiCore'
	
    conda params.soloTE_conda

    input:
        tuple val(name), path(bam)
        val (outputdir)
        val (te_bed)
        val (solote_dir)

    output:
        tuple val(name), path ("${name}_SoloTE_output"), emit: matrix

    publishDir "${outputdir}",
        mode: "link", overwrite: true, enabled: !params.no_output

    script:
        cores = 8
        """
        python ${solote_dir}/SoloTE_pipeline.py \\
            --threads ${cores} \\
            --bam ${bam} \\
            --teannotation ${te_bed} \\
	    --outputdir . \\
            --outputprefix ${name}
        """
}
