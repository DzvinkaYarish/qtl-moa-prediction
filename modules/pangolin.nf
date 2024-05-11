process GetPangolinPreds {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/pt/"
    label 'gpu'

    input:
    path infile
    val cell_type

    output:
    path 'pangolin_preds.vcf'

    shell:
    '''
    export workdir=$(pwd)

    module load any/python/3.8.3-conda

    conda activate pt

    !{params.HOME}/.conda/envs/pt/bin/python -c "import torch; print(torch.cuda.is_available())"

    !{params.HOME}/.conda/envs/pt/bin/python !{params.HOME}/Thesis/pipelines/python_scripts/convert_to_pangolin_format.py !{infile} pangolin_input.vcf

    !{params.HOME}/.conda/envs/pt/bin/pangolin ${workdir}/pangolin_input.vcf \
               /gpfs/space/projects/genomic_references/annotations/hg38/hg38.fa \
                !{params.HOME}/Thesis/common_data/gencode.v44.basic.annotation.db \
                pangolin_preds \
                -d 500
    '''
 }