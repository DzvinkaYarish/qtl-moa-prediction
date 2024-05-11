process build_hisat2_index{
    container = 'quay.io/eqtlcatalogue/rnaseq_hisat2:v22.03.01'

    input:
    file(reference_genome)

    output:
    tuple file("hisat2_index.*.ht2")

    script:
    """
    hisat2-build ${reference_genome} hisat2_index
    """
}


process hisat2_align{
    container = 'quay.io/eqtlcatalogue/rnaseq_hisat2:v22.03.01'
    publishDir "/gpfs/space/home/dzvenymy/Bioinformatics/hw3/results/"

    input:
    tuple val(sample_name), file(fastq1), file(fastq2)
    file hisat2_indices

    output:
    tuple val(sample_name), file("${sample_name}.sortedByCoords.bam"), file("${sample_name}.sortedByCoords.bam.bai")

    script:
    index_base = hisat2_indices[0].toString() - ~/.\d.ht2/
    """
    hisat2 -x $index_base -1 ${fastq1} -2 ${fastq2} | samtools view -Sb > ${sample_name}.bam
    samtools sort -o ${sample_name}.sortedByCoords.bam ${sample_name}.bam
    samtools index ${sample_name}.sortedByCoords.bam
    """
}

// process samsort{
//
//  input:
//  tuple val(sample_name), file(bam)
//
//  output:
//  tuple val(sample_name), file("${sample_name}.sortedByCoords.bam")
//
//  script:
//  """
//
//  """
// }

process feature_count{
    container = 'docker://quay.io/eqtlcatalogue/rnaseq:v20.11.1'
    publishDir "/gpfs/space/home/dzvenymy/Bioinformatics/hw3/results/"
    input:
    file gtf
    tuple val(sample_name), file(bam), file(bai)

    output:
    file("${sample_name}.counts")

    script:
    """
    featureCounts -p -C -D 5000 -d 50 -s2 -a ${gtf} -o ${sample_name}.counts ${bam}
    """
}