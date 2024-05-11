process FindOverlaps {
    container = 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3'
    input:
    path bed
    path infile


    output:
    path "${infile.simpleName}_${bed.simpleName}_overlap.bed"

    shell:
    '''
    awk -F'\\t' 'NR>1 {

    chromosome=a[1];
    position=a[2];
    start=position;
    end=position;
    print chromosome"\\t"start"\\t"end
    }' !{infile} > !{infile.simpleName}.bed

    awk '!seen[$0]++' !{infile.simpleName}.bed > !{infile.simpleName}_unique.bed

    bedtools intersect -a !{infile.simpleName}_unique.bed -b !{bed} > !{infile.simpleName}_!{bed.simpleName}_overlap.bed
    '''
}

process CombineOverlaps {
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    path overlaps
    val outfile_name

    output:
    path "${outfile_name}.bed"

    script:
    """
    cat ${overlaps.join(' ')} > ${outfile_name}.bed
    """
}