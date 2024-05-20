This repo contains the code for the  Master's thesis 'Predicting the molecular mechanisms of genetic variants'.

It implements a Nextflow workflow, which retrieves all necessary genomic features and runs inference with several neural networks to obtain features for the mode of action prediction model - classifying whether the (fine-mapped) variant acts as a spicing (s)QTL or gene expression affected by chromatin accessibility (ce)QTL.

#### Genomic features
1. Binary variable indicating whether the variant is located within the gene body
2. Distance from the variant to the closest annotated splice junction (GENCODE v39
annotation)
3. Number of overlaps with open chromatin regions in 5 cell types. We used the same
ENCODE DNASE/ATAC-seq experiments on which ChromBPNet models were
trained.
4. Number of overlaps with binding sites of RNA binding proteins. We took the
binding sites of 211 RBPs, identified by Nostrand et al.

#### Neural features
1. Splicing scores from SpliceAI and Pangolin. Each model produces two scores:
maximum increase and decrease in the probability of a site being a splice junction
in a 1000bp window around the variant.
2. Enformer SAD scores for five CAGE tracks (gene expression) and five DNASE
tracks. SAD score is a difference between Enformer predictions for reference
and alternative alleles, averaged over the eight flanking bins representing 1000bp
window.
3. ChromBPNet difference scores for five cell types.
