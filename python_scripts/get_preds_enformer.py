import h5py
import pandas as pd
import numpy as np
from liftover import get_lifter
import sys

chr1 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.1.h5"
chr2 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.2.h5"
chr3 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.3.h5"
chr4 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.4.h5"
chr5 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.5.h5"
chr6 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.6.h5"
chr7 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.7.h5"
chr8 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.8.h5"
chr9 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.9.h5"
chr10 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.10.h5"
chr11 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.11.h5"
chr12 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.12.h5"
chr13 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.13.h5"
chr14 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.14.h5"
chr15 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.15.h5"
chr16 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.16.h5"
chr17 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.17.h5"
chr18 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.18.h5"
chr19 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.19.h5"
chr20 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.20.h5"
chr21 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.21.h5"
chr22 = "/gpfs/space/projects/genomic_references/enformer/1000G.MAF_threshold=0.005.22.h5"

converter = get_lifter('hg38', 'hg19')


def lift(row):
    try:
        return converter[row['CHR']][row['BP_hg38']][0][1]
    except IndexError:
        # print(converter[row['CHR']][row['BP_hg38']])
        return -1


def get_chr(x):
    ch = x.split('_')[0].strip('chr')
    try:
        return int(ch)
    except ValueError:
        print(x)
        return ch


def format_input(df):
    f_df = pd.DataFrame()
    f_df['SNP'] = df['rsid']
    f_df['A1'] = df['variant'].apply(lambda x: x.split('_')[2])
    f_df['A2'] = df['variant'].apply(lambda x: x.split('_')[3])
    f_df['BP_hg38'] = df['variant'].apply(lambda x: x.split('_')[1]).astype(int)
    f_df['CHR'] = df['variant'].apply(get_chr)
    f_df['BP_hg19'] = f_df[['BP_hg38', 'CHR']].apply(lift, axis=1)
    f_df['variant'] = df['variant']

    return f_df


def compare(f, comp_chr, enf_pos, enf_snp, cell_track_id):
    comparison_df = pd.DataFrame(
        columns=['rsid', 'variant', 'SAD_selected', 'ref_orig', 'alt_orig', 'ref_Enformer', 'alt_Enformer', 'pos_orig',
                 'pos_Enformer'])

    for i, row in comp_chr.iterrows():

        rsid = row['SNP']
        pos = row['BP_hg19']
        match_snp = np.where(enf_snp == rsid)[0]
        match_pos = np.where(enf_pos == pos)[0]
        if len(match_snp) + len(match_pos) > 0:
            indx = match_snp[0] if len(match_snp) > 0 else match_pos[0]
            sad = f['SAD'][indx][cell_track_id]
            ref = f['ref'][indx].astype(str)
            alt = f['alt'][indx].astype(str)
            p = f['pos'][indx].astype(str)
        else:
            sad = -1000
            ref = row['A1']
            alt = row['A2']
            p = row['BP_hg19']

        comparison_df.loc[i] = [rsid,
                                row['variant'],
                                sad,
                                row['A1'],
                                row['A2'],
                                ref,
                                alt,
                                row['BP_hg19'],
                                p
                                ]
    return comparison_df


def check_ref_alt(df):
    alt_ref_switched = np.array(df.loc[(df["ref_orig"].values == df["alt_Enformer"].values) & (
                df["alt_orig"].values == df["ref_Enformer"].values)].index)
    # print(alt_ref_switched)
    df.loc[alt_ref_switched, 'SAD_selected'] = df.loc[alt_ref_switched, 'SAD_selected'] * (-1)


def get_enformer_preds(comp, cell_type_indx):

    Enformer_predictions = [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14,
                            chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22]

    total_dfs = pd.DataFrame(
        columns=['rsid', 'variant', 'SAD_selected', 'ref_orig', 'alt_orig', 'ref_Enformer', 'alt_Enformer', 'pos_orig',
                 'pos_Enformer'])

    c = 1
    for file in Enformer_predictions:

        comp_chr = comp[comp.CHR == c]

        if len(comp_chr) > 0:
            f = h5py.File(file, "r")
            pos = np.array(f['pos'])
            snp = np.array(f['snp']).astype(str)

            comp_df = compare(f, comp_chr, pos, snp, cell_type_indx)

            check_ref_alt(comp_df)

            total_dfs = pd.concat([total_dfs, comp_df], ignore_index=True)

            c += 1
    total_dfs = total_dfs[['rsid', 'variant', 'SAD_selected']]
    return total_dfs


if __name__ == '__main__':
    input_file = sys.argv[1]
    input_df = pd.read_csv(input_file, sep='\t')

    cell_type = int(sys.argv[3]) # ToDo: Add cell type resolution

    input_df = format_input(input_df)
    out_df = get_enformer_preds(input_df, cell_type)

    out_file = sys.argv[2]
    out_df.to_csv(out_file, sep='\t', index=False)
