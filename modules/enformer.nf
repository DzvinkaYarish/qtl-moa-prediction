process GetEnformerPreds {
    publishDir "${params.out_dir}", mode: 'copy'
    conda "${params.HOME}/.conda/envs/chrombpnet/"
    label 'gpu'

    input:
    path infile
    val cell_type

    output:
    path "${infile.simpleName}_enformer_preds.csv"

    script:
    """

    module purge
    module load cuda/11.7.0
    module load cudnn/8.2.0.53-11.3
    module load any/python/3.8.3-conda

    conda activate chrombpnet

    ${params.HOME}/.conda/envs/chrombpnet/bin/python ${params.HOME}/Thesis/pipelines/python_scripts/get_enformer_preds_via_inference.py  $infile ${infile.simpleName}_enformer_preds.csv $cell_type
    """
}