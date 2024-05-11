process GetGenomicDistances {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/qtl/"
    input:
    path 'infile.csv'
    val cell_type

    output:
    path 'dist_features.csv'

    script:
    """
    module load any/python/3.8.3-conda
    conda activate qtl
    ${params.HOME}/.conda/envs/qtl/bin/python ${params.HOME}/Thesis/pipelines/python_scripts/get_preds_enformer.py  infile.csv enformer_preds.csv $cell_type
    """
}