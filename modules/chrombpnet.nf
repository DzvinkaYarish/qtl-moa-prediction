process GetChrombpnetPreds {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/chrombpnet/"
    label 'gpu'

    input:
    path infile
    val cell_type

    output:
    path "chrompbnet_preds_${cell_type}.csv"

    stub:
    """
    cat /gpfs/space/home/dzvenymy/chrombpnet_tutorial/variant_scoring/Bossini-Castillo_2019_QTD100008/ENCSR000EPK_model/variant_scores.tsv > variant_scores.tsv
    """

    shell:
    '''
    export workdir=$(pwd)
    echo $workdir
    module purge
    module load cuda/11.7.0
    module load cudnn/8.2.0.53-11.3
    module load any/python/3.8.3-conda

    conda activate chrombpnet

    !{params.HOME}/.conda/envs/chrombpnet/bin/python -c  "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

    !{params.HOME}/.conda/envs/chrombpnet/bin/python !{params.HOME}/Thesis/pipelines/python_scripts/convert_to_chrombpnet_format.py !{infile} 'chrombpnet_preds.csv'

    !{params.HOME}/.conda/envs/chrombpnet/bin/python !{params.HOME}/variant-scorer/src/variant_scoring.py \
    -l ${workdir}/chrombpnet_preds.csv \
    -g !{params.HOME}/chrombpnet_tutorial/common_data/hg38.fa \
    -m !{params.HOME}/chrombpnet_tutorial/!{cell_type}/chrombpnet_model/models/chrombpnet_nobias.h5 \
    -o ${workdir}/chrompbnet_preds_!{cell_type}.csv \
    -s !{params.HOME}/chrombpnet_tutorial/common_data/hg38.chrom.sizes \
    --no_hdf5 \
    -n 1
    '''

}