nextflow.enable.dsl=2

params.no_output = false


//script_path="/bi/group/bioinf/Laura_B/nf_yang_wang/scripts/"
//splitpool_path = "/bi/group/bioinf/Laura_B/nf_yang_wang/splitpoolquantitation/"

// ONLY KEEP UNIQUE READS - where is this used? - just in the main workflow
process SAMTOOLS_FILT {

    // 	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
        tuple val(name), path (bam)
		val (outputdir)
		val (samtools_filt_args)	

	output:
		tuple val(name), path ("*bam"), 	  emit: bam

    publishDir "${outputdir}",
		mode: "link", overwrite: true
		
		"""
		module load samtools || echo "no module found"
		samtools view -h -q 10 $bam -o ${bam}_clean.bam
        rename .bam_clean _clean *
		"""	
}



// process SAMTOOLS_FILT {

//     // 	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
// 	label 'bigMem' // 20GB

// 	input:
//         tuple val(name), path (bam)
// 		//path(bam)	
//         val (outputdir)	

// 	output:
// 		tuple val(name), path ("*bam"), 	  emit: bam

//     publishDir "${outputdir}/nf_chosen_outputs",
// 		mode: "link", overwrite: true
		
// 		"""
// 		module load samtools || echo "no module found"
// 		samtools view -h -q 10 $bam | samtools sort -o ${bam}_clean_sorted.bam
//         rename .bam_clean _clean *
// 		"""	
// }





process SAMTOOLS_SORT{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		tuple val(name), path (bam)
		val (outputdir)
		val (samtools_sort_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		tuple val(name), path ("*bam"), 	  emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	
    script:
		samtools_sort_options = samtools_sort_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}
		
		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load samtools || echo "no module found"
		samtools sort $samtools_sort_options $bam -o ${bam}_sorted.bam 
		rename .bam_sorted _sorted *
    	"""	
}


// output is the index file and the sorted bam so that they stay together
process SAMTOOLS_INDEX{	
    
	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		tuple val(name), path (bam)
		val (outputdir)
		val (samtools_index_args)
		val (verbose)

	output:
		path "*.bai",     emit: bai
		tuple val(name), path (bam), 	  emit: bam
    	
	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

    script:
		samtools_index_options = samtools_index_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS INDEX ARGS: " + samtools_index_args)
		}
		
		"""
		module load samtools || echo "no module found"
		samtools index $samtools_index_options $bam
		"""	
}

process SAMTOOLS_SORT_SAM{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		tuple val(name), path (sam)
		val (outputdir)
		val (samtools_sort_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		tuple val(name), path ("*bam"), 	  emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	
    script:
		samtools_sort_options = samtools_sort_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}
		
		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load samtools || echo "no module found"
		samtools view $sam | samtools sort $samtools_sort_options -o ${sam}_sorted.bam 
		rename .bam_sorted _sorted *
    	"""	
}