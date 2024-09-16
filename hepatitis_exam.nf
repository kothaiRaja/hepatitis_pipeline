nextflow.enable.dsl = 2

params.accession = "M21012"
params.store = "${launchDir}/cache"
params.out = "$launchDir/output"
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
		path infile

    output:
		path "*"

    script:
    """
    cat ${infile} > combined.fasta
    """
	
}

workflow{
	downloadchannel = download_reference(Channel.from(params.accession))
	input_fasta_channel = Channel.fromPath("${params.in}/*.fasta")
	combined_channel = downloadchannel.concat(input_fasta_channel)
	combined_fastafiles = combine_fasta(combined_channel.collect())
}