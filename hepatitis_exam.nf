nextflow.enable.dsl = 2

//params.accession: This stores the NCBI accession number (M21012) to be used to download a reference sequence.
params.accession = "M21012"
//The directory where intermediate files (such as the downloaded reference sequence) are stored. ${launchDir} refers to the directory where the Nextflow script is run.
params.store = "${launchDir}/cache"
//The output directory where the final results are published.
params.out = "$launchDir/output"
//The directory where the input FASTA files (from local sequencing results) are located.
params.in = "${launchDir}/hepatitis_fasta"
//downloading reference sequence

process download_reference{
	storeDir params.store
	input:
		val accession
	output:
		path "M21012.fasta"
	"""
	wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M21012&rettype=fasta&retmode=text" -O M21012.fasta
	"""
	
}


process combine_fasta {
    publishDir params.out, mode: "copy", overwrite: true

    input:
		path infile //The process takes a list of FASTA files (both the downloaded reference and the local sequencing files

    output:
		path "combined.fasta"

    script:
    """
    cat ${infile} > combined.fasta
    """	
}

//Multiple Alignment using Fast fourier transform

process Aligner_MAFFT {
	publishDir params.out, mode: "copy", overwrite: true
	container "https://depot.galaxyproject.org/singularity/mafft%3A7.525--h031d066_1"
	input:
		path infile
	output:
		path "mafft.fasta"
	"""
	mafft --auto ${infile} > mafft.fasta 
	"""
}

process trimAll{
	publishDir params.out, mode: "copy", overwrite: true
	container "https://depot.galaxyproject.org/singularity/trimal%3A1.5.0--h4ac6f70_1"
	input:
		path infile
	output: 
		path "*"
		
//trimal is run with the -automated1 option to automatically clean up the alignment and generate both the trimmed FASTA and the HTML report.
	"""
	trimal -in ${infile} -out alignment_ref-seq_${params.accession}.fasta -htmlout report_ref-seq_${params.accession}.html -automated1 
	"""
}
workflow{
	downloadchannel = download_reference(Channel.from(params.accession))
	input_fasta_channel = Channel.fromPath("${params.in}/*.fasta")
	combined_channel = downloadchannel.concat(input_fasta_channel)
	combined_fastafiles = combine_fasta(combined_channel.collect())
	Alignerchannel = Aligner_MAFFT(combined_fastafiles)
	trimchannel = trimAll(Alignerchannel)
}