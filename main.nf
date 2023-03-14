/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
params.index = "index"


log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    
    publishDir params.index, mode: 'copy'

    input:
    path transcriptome
    //we are going to pass a channel to Index; and each file from that will be an input variable called 
    //transcriptome.  In this case, there is only one file.
    output:
    path 'salmon_index'
    //the output will be a file that we assign to the variable 'salmon_index'.
    //in this case, we know the output will be just one file
    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {

    tag "$sample_id" //for traceability in the log
    publishDir params.outdir, mode: 'copy' //for writing data into our results folder

    input:
    path salmon_index 
    //first channel input is file coming from index_ch channel.  There will only be one file.
    tuple val(sample_id), path(reads)
    //second channel input is a tuple that contains the sample_id (val) and then the two files 

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    //create a channel that retrieves a tuple where the first element is the read pair prefix,
    //and the second element is a list of the full path file names

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)   
}
