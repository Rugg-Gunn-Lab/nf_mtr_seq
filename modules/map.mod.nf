nextflow.enable.dsl=2
params.local = ''
params.no_output = false

process BOWTIE2 {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)

	output:
	    tuple val(name), path ("*bam"),        emit: bam
		path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		cores = 8
		readString = ""

		// Options we add are
		bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file
		
		if (params.local == '--local'){
			// println ("Adding option: " + params.local )
			bowtie_options += " ${params.local} " 
		}

		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}


		index = params.genome["bowtie2"]
		bowtie_name = name + "_" + params.genome["name"]

		"""
		module load bowtie2 || echo "no module found"
		module load samtools || echo "no module found"
		bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}



process STAR {

    label 'hugeMem'
    label 'multiCore'

    input:
        tuple val(name), path(reads)
        val (outputdir)
        val (star_index)

	output:
	    tuple val(name), path ("*bam"),  emit: bam
		//path "*final.out",               emit: stats 
        path "*.txt",               emit: stats 

    //publishDir "${outputdir}/nf_chosen_outputs",
	publishDir "${outputdir}",
		mode: "link", overwrite: true

    script:
       
        """
        module load star || echo "no module found"
        module load samtools || echo "no module found"
        STAR --runThreadN 16 --genomeDir ${star_index} --readFilesIn ${reads} --readFilesCommand zcat --outStd SAM 2>${name}_star_out2.txt | samtools view -h -b -q 50 --threads 16 - | samtools sort --threads 16 -o ${name}_filt_sorted.bam
        """
}

// for mapping to the barcode genome
process BOWTIE_BC {

    label 'bigMem' // 20GB
    label 'multiCore'

    input:
        tuple val(name), path(reads)
        val (outputdir)
        val (bowtie_args)
        val (barcode_genome)

	output:
	    tuple val(name), path ("*sam"),     emit: sam
		path "*stats_barcodes.txt",         emit: stats 

    // we only want the stats file written out as a hard link, not the sam file
    publishDir = [
        //path: { "${outputdir}/nf_chosen_outputs" },
		path: { "${outputdir}" },
        mode: "link", overwrite: true,
        pattern: "*.txt"
    ]

    script:
        cores = 8

        """
        module load bowtie || echo "no module found"
        module load samtools || echo "no module found"
        bowtie ${barcode_genome} --sam -p ${cores} ${bowtie_args} ${reads} 2>${name}_bowtie_stats_barcodes.txt ${name}_mapped_BC.sam
        """
}