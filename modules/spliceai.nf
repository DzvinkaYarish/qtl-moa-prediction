process GetSpliceAIPreds {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/tf1/"
    label 'gpu'

    input:
    path infile
    val cell_type

    output:
    path 'spliceai_preds.vcf'

    shell:
    '''
    pwd
    export workdir=$(pwd)
    echo $workdir


    module load cuda/10.0.130 && module load cudnn/7.6.5.32-10.0

    module load any/python/3.8.3-conda

    conda activate tf1


    !{params.HOME}/.conda/envs/tf1/bin/python -c  "import tensorflow as tf; print(tf.test.is_gpu_available())"




    !{params.HOME}/.conda/envs/tf1/bin/python !{params.HOME}/Thesis/pipelines/python_scripts/convert_to_pangolin_format.py !{infile} spliceai_input.vcf

    !{params.HOME}/.conda/envs/tf1/bin/spliceai -I ${workdir}/spliceai_input.vcf \
                                                -O spliceai_preds.vcf \
                                                -R /gpfs/space/projects/genomic_references/annotations/hg38/hg38.fa \
                                                -A !{params.HOME}/Thesis/spliceai_annotations/gencode.v44.basic.annotation.txt.gz \
                                                -D 500 2>&1
    '''
}