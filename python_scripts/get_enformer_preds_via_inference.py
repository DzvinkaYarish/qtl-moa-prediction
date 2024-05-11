import tensorflow as tf
import tensorflow_hub as hub

import tensorflow_hub as hub
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
from tqdm import tqdm

import sys

root = '/gpfs/space/home/dzvenymy'
transform_path = root + '/qtl_labeling/enformer_data/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file = '/gpfs/space/projects/genomic_references/annotations/hg38/hg38.fa'


# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)


class EnformerScoreVariantsRaw:

  def __init__(self, tfhub_url, organism='human'):
    self._model = Enformer(tfhub_url)
    self._organism = organism

  def predict_on_batch(self, inputs, crop=False, **kwargs):
    ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
    alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]

    out_len = ref_prediction.shape[1]
    if crop:
        print('Crop')
        alt_prediction = alt_prediction[:, int(out_len / 2) - 4: int(out_len / 2) + 4]
        ref_prediction = ref_prediction[:, int(out_len / 2) - 4: int(out_len / 2) + 4]

    return alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)


class EnformerScoreVariantsNormalized:

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      transform_pipeline = joblib.load(f)
    self._transform = transform_pipeline.steps[0][1]  # StandardScaler.

  def predict_on_batch(self, inputs, **kwargs):
    scores = self._model.predict_on_batch(inputs, **kwargs)
    return self._transform.transform(scores)


class EnformerScoreVariantsPCANormalized:

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      self._transform = joblib.load(f)
    self._num_top_features = num_top_features

  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)[:, :self._num_top_features]

# @title `variant_centered_sequences`

class FastaStringExtractor:

    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


def variant_generator(vcf_file, gzipped=False):
  """Yields a kipoiseq.dataclasses.Variant for each row in VCF file."""
  def _open(file):
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file)

  with _open(vcf_file) as f:
    for line in f:
      if line.startswith('#'):
        continue
      chrom, pos, id, ref, alt_list = line.split('\t')[:5]
      # Split ALT alleles and return individual variants as output.
      for alt in alt_list.split(','):
        yield kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos,
                                           ref=ref, alt=alt, id=id)


def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)


def variant_centered_sequences(vcf_file, sequence_length, gzipped=False,
                               chr_prefix=''):
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))

  for variant in variant_generator(vcf_file, gzipped=gzipped):
    interval = Interval(chr_prefix + variant.chrom,
                        variant.pos, variant.pos)
    interval = interval.resize(sequence_length)
    center = interval.center() - interval.start

    reference = seq_extractor.extract(interval, [], anchor=center)
    alternate = seq_extractor.extract(interval, [variant], anchor=center)

    yield {'inputs': {'ref': one_hot_encode(reference),
                      'alt': one_hot_encode(alternate)},
           'metadata': {'chrom': chr_prefix + variant.chrom,
                        'pos': variant.pos,
                        'id': variant.id,
                        'ref': variant.ref,
                        'alt': variant.alt}}


if __name__ == '__main__':
    model = EnformerScoreVariantsNormalized(model_path, transform_path)
    fasta_extractor = FastaStringExtractor(fasta_file)
    input_file = sys.argv[1]
    var_df = pd.read_csv(input_file, sep='\t')

    ref_data = []
    alt_data = []
    invalid_vars = []
    for i, row in tqdm(var_df.iterrows()):
        ch, pos, ref, alt = row['variant'].split('_')
        variant = kipoiseq.Variant(ch, pos, ref, alt)

        # Center the interval at the variant
        interval = kipoiseq.Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
        seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
        center = interval.center() - interval.start

        try:
            reference = seq_extractor.extract(interval, [], anchor=center)
            alternate = seq_extractor.extract(interval, [variant], anchor=center)

            ref_data.append(one_hot_encode(reference))
            alt_data.append(one_hot_encode(alternate))
        except KeyError:
            print(row['variant'])
            invalid_vars.append(i)
            continue
        except ValueError as exp:
            print(exp)
            print(row['variant'])
            invalid_vars.append(i)

            continue

    ref_data = np.stack(ref_data)
    alt_data = np.stack(alt_data)

    bs = 10
    cage_tracks = [5110, 5213, 5218, 4758]
    dnase_tracks = [14, 62, 156, 166, 41]
    tracks = cage_tracks + dnase_tracks
    columns = [f'SAD_score_{i}' for i in tracks]
    sad_scores = pd.DataFrame(columns=columns)
    for i in tqdm(range(0, len(ref_data), bs)):
        if i + bs > len(ref_data):
            ref_batch = ref_data[i:]
            alt_batch = alt_data[i:]
        else:
            ref_batch = ref_data[i:i + bs]
            alt_batch = alt_data[i:i + bs]
        preds = model.predict_on_batch({'ref': ref_batch, 'alt': alt_batch}, crop=True)
        sad_scores = pd.concat((sad_scores, pd.DataFrame(data=preds[:, tracks], columns=columns)))
    var_df = var_df.drop(index=invalid_vars).reset_index(drop=True)
    var_df = pd.concat((var_df, sad_scores.reset_index(drop=True)), axis=1)
    out_file = sys.argv[2]
    var_df.to_csv(out_file, sep='\t', index=False)

