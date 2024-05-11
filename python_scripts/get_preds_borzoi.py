import json
import sys
from tqdm import tqdm

import numpy as np
import pandas as pd
import pysam
import pyfaidx
import tensorflow as tf

from baskerville import seqnn
from baskerville import gene as bgene
from baskerville import dna

from borzoi_helpers import *

tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

pyfaidx.Faidx('/gpfs/space/projects/genomic_references/annotations/hg38/hg38.fa')

# Model configuration

params_file = '/gpfs/space/home/dzvenymy/borzoi/examples/params_pred.json'
targets_file = 'targets_gtex.txt'  # Subset of targets_human.txt

seq_len = 524288
n_folds = 1  # To use only one model fold, set to 'n_folds = 1'. To use all four folds, set 'n_folds = 4'.
rc = False  # Average across reverse-complement prediction
tracks_ids = [6112,
 6113,
 6114,
 6115,
 6116,
 6117,
 6118,
 6119,
 6120,
 6121,
 6204,
 6205,
 6340,
 6434,
 6514,
 6515,
 6668,
 6669,
 7306,
 7307,
 7334,
 7394,
 1288,
 6517,
 6518,
 6864,
 6865,
 7460,
 7461,
 7514,
 7515,
 1338,
 6574,
 7106,
 7297,
 1432,
 6757,
 6758,
 7219,
 7454,
 7455,
 1442,
 6292,
 6293,
 6765,
 6873,
 7423,
 7504,
 1317]

# Read model parameters

with open(params_file) as params_open:
    params = json.load(params_open)

    params_model = params['model']
    params_train = params['train']

# Read targets

# targets_df = pd.read_csv(targets_file, index_col=0, sep='\t')
# target_index = targets_df.index

target_index = tracks_ids

# Create local index of strand_pair (relative to sliced targets)
if rc:
    strand_pair = targets_df.strand_pair

    target_slice_dict = {ix: i for i, ix in enumerate(target_index.values.tolist())}
    slice_pair = np.array([
        target_slice_dict[ix] if ix in target_slice_dict else ix for ix in strand_pair.values.tolist()
    ], dtype='int32')

# Initialize model ensemble

models = []
for fold_ix in range(n_folds):

    model_file = "/gpfs/space/home/dzvenymy/borzoi/examples/saved_models/f" + str(fold_ix) + "/model0_best.h5"

    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(model_file, 0)
    seqnn_model.build_slice(target_index)

    if rc:
        seqnn_model.strand_pair.append(slice_pair)

    seqnn_model.build_ensemble(rc, [0])

    models.append(seqnn_model)

fasta_open = pysam.Fastafile('/gpfs/space/projects/genomic_references/annotations/hg38/hg38.fa')

transcriptome = bgene.Transcriptome('/gpfs/space/home/dzvenymy/Thesis/common_data/gencode.v39.annotation.nochr.gtf')


if __name__ == '__main__':
    input_file = sys.argv[1]
    var_df = pd.read_csv(input_file, sep='\t')

    preds_df = pd.DataFrame()

    for i, row in tqdm(var_df.iterrows()):
        chrom, pos, ref, alt = row['variant'].split('_')
        # if i > 1:
        #     if len(preds_df[preds_df['variant'] == row['variant']]) > 0 :
        #         preds_df = pd.concat([preds_df, preds_df[preds_df['variant'] == row['variant']].iloc[0]])
        #         continue
        pos = int(pos)
        start = pos - seq_len // 2
        end = pos + seq_len // 2

        search_gene = row['gene_id']
        if len(alt) == 1 and isinstance(search_gene, str): # awkward way to filter out genes which are nones
            print(f'processing {row["variant"]}...')
            gene_keys = [gene_key for gene_key in transcriptome.genes.keys() if search_gene in gene_key]

            gene = transcriptome.genes[gene_keys[0]]

            # Determine output sequence star
            seq_out_start = start + seqnn_model.model_strides[0] * seqnn_model.target_crops[0]
            seq_out_len = seqnn_model.model_strides[0] * seqnn_model.target_lengths[0]

            # Determine output positions of gene exons
            gene_slice = gene.output_slice(seq_out_start, seq_out_len, seqnn_model.model_strides[0], True)

            poses = [pos]
            alts = [alt]

            sequence_one_hot_wt = process_sequence(fasta_open, chrom, start, end)

            #Induce mutation(s)
            sequence_one_hot_mut = np.copy(sequence_one_hot_wt)

            for pos, alt in zip(poses, alts) :
                alt_ix = -1
                if alt == 'A' :
                    alt_ix = 0
                elif alt == 'C' :
                    alt_ix = 1
                elif alt == 'G' :
                    alt_ix = 2
                elif alt == 'T' :
                    alt_ix = 3

                sequence_one_hot_mut[pos-start-1] = 0.
                sequence_one_hot_mut[pos-start-1, alt_ix] = 1.

            #Make predictions
            y_wt = predict_tracks(models, sequence_one_hot_wt)
            y_mut = predict_tracks(models, sequence_one_hot_mut)

            y_mut = y_mut[0][0]
            y_wt = y_wt[0][0]

            ## compute scores

            log_scores = np.sqrt(np.sum(np.power(np.log2(1 + y_mut) - np.log2(1 + y_wt), 2), axis=0))
            sum_scores = np.mean(y_mut - y_wt, axis=0)

            splice_scores = np.max(np.abs(y_mut[gene_slice] / np.sum(y_mut[gene_slice], axis=0) - y_wt[gene_slice] / np.sum(y_wt[gene_slice], axis=0)),axis=0)
        else:
            print('skipping variant...')
            log_scores = [0]*len(tracks_ids)
            sum_scores = [0]*len(tracks_ids)
            splice_scores = [0]*len(tracks_ids)

        res_df = pd.concat((pd.DataFrame({'variant': [row['variant']]}),
                            pd.DataFrame({f'Borzoi_{i}_log_score': [v] for i,v in zip(tracks_ids, log_scores)}),
                           pd.DataFrame({f'Borzoi_{i}_sum_score': [v] for i,v in zip(tracks_ids, sum_scores)}),
                        pd.DataFrame({f'Borzoi_{i}_splice_score':[v] for i,v in zip(tracks_ids, splice_scores)}),

                           ), axis=1)
        preds_df = pd.concat((preds_df, res_df), axis=0)

    out_file = sys.argv[2]
    preds_df.to_csv(out_file, sep='\t', index=False)
