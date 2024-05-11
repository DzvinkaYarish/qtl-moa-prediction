process GetFeatures {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/qtl/"
    input:
    path all_preds

    output:
    path 'features.csv'

    script:
    """
    module load any/python/3.8.3-conda
    conda activate qtl
    ${params.HOME}/.conda/envs/qtl/bin/python ${params.HOME}/Thesis/pipelines/python_scripts/prepare_features.py ${all_preds.join(' ')}
    """
}