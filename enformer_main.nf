nextflow.enable.dsl=2

include { GetEnformerPreds } from './modules/enformer'

workflow {
    infile = "${params.in_file}"
    cell_type = "${params.cell_type}"
//     datasets = Channel.fromPath( "/gpfs/space/home/dzvenymy/qtl_labeling/cQTL_data/AFGR*.csv" )

    out1 = GetEnformerPreds(infile, cell_type)

}
