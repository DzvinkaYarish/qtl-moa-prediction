import pandas as pd
import sys

import warnings
warnings.filterwarnings('ignore')


def format_spliceai(fpath, cols_to_use=['spliceai_increase', 'spliceai_decrease']):
    cols = [str(i) for i in range(26)]
    cols[0] = "chr"
    cols[1] = "pos"

    preds_df = pd.read_csv(fpath,
        sep=r'\t|\|\s*', header=None, comment='#', engine='python', names=cols,
        on_bad_lines=lambda bad_line: bad_line[:26])
    preds_df = preds_df.drop(["2", "3", "4", "5", "6", "7", "8"], axis=1)

    preds_df = preds_df.rename({v: k for k, v in
                                zip(['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL'],
                                    map(str, list(range(9, 17))))}, axis=1)
    preds_df.loc[preds_df['DS_AG'] == '.', ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL' ]] = [0, 0, 0, 0, 0, 0, 0, 0]
    preds_df[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']] = preds_df[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].astype('float32')

    preds_df['spliceai_increase'] = preds_df["DS_AG"] + preds_df["DS_DG"]
    preds_df['spliceai_decrease'] = preds_df["DS_AL"] + preds_df["DS_DL"]
    preds_df = preds_df.reset_index(drop=True)

    return preds_df[cols_to_use]


def format_pangolin(fpath, cols_to_use=['pangolin_increase', 'pangolin_decrease']):
    cols = [str(i) for i in range(13)]
    cols[0] = "chr"
    cols[1] = "pos"

    preds_df = pd.read_csv(fpath, sep="\t|\||:", header=None, comment='#',
                           engine='python', names=cols,
                           on_bad_lines=lambda bad_line: bad_line[:13], index_col=False)
    preds_df = preds_df.drop(["2", "3", "4", "5", "6", "7", "12"], axis=1)

    preds_df = preds_df.rename(
        {v: k for k, v in zip(["pos_increase", "pangolin_increase", "pos_decrease", "pangolin_decrease"], map(str, list(range(8, 12))))},
        axis=1)
    preds_df = preds_df.reset_index(drop=True)

    return preds_df[cols_to_use]


def format_chrombpnet(fpath, prefix='', cols_to_use=['logfc']):
    preds_df = pd.read_csv(fpath, sep='\t')
    preds_df = preds_df.rename(columns={col: f'{prefix}{col}' for col in cols_to_use})
    return preds_df[[f'{prefix}{col}' for col in cols_to_use]]


def format_rbp(fpath, df):
    rbps = pd.read_csv(fpath, sep='\t', names=['CHR', 'START', 'END'])
    rbps['pos'] = rbps.apply(lambda x: x['CHR'] + '_' + str(x['START']), axis=1)
    rbps['# RBP sites'] = rbps['pos'].map(rbps['pos'].value_counts())
    rbps = rbps.drop_duplicates()
    df['pos'] = df['variant'].apply(lambda x: '_'.join(x.split('_')[:2]))
    df = pd.merge(df, rbps[['pos', '# RBP sites']], how='left', on='pos')
    return df

def format_peaks(fpath, df):
    peaks = pd.read_csv(fpath, sep='\t', names=['CHR', 'START', 'END'])
    peaks['pos'] = peaks.apply(lambda x: x['CHR'] + '_' + str(x['START']), axis=1)
    peaks['# peak overlaps'] = peaks['pos'].map(peaks['pos'].value_counts())
    peaks = peaks.drop_duplicates()
    df['pos'] = df['variant'].apply(lambda x: '_'.join(x.split('_')[:2]))
    df = pd.merge(df, peaks[['pos', '# peak overlaps']], how='left', on='pos')
    return df


if __name__ == '__main__':
    # enformer
    # spliceai
    # pangolin
    # rbp
    # peaks
    # chrombpnet

    enf_preds = pd.read_csv(sys.argv[1], sep='\t')[['variant',  'SAD_score_5110', 'SAD_score_5213', 'SAD_score_5218', 'SAD_score_4758', 'SAD_score_14', 'SAD_score_62', 'SAD_score_156', 'SAD_score_166', 'SAD_score_41']]

    spliceai_preds = format_spliceai(sys.argv[2])
    pangolin_preds = format_pangolin(sys.argv[3])
    features_df = pd.concat((enf_preds, spliceai_preds, pangolin_preds), axis=1)

    features_df = format_rbp(sys.argv[4], features_df)
    features_df = format_peaks(sys.argv[5], features_df)

    for fpath in sys.argv[6:]:
        chrombpnet_preds = format_chrombpnet(fpath, prefix=fpath.split('/')[-1].split('.')[0].split('_')[-1] + '_')
        features_df = pd.concat((features_df, chrombpnet_preds), axis=1)

    features_df.to_csv('features.csv', sep='\t', index=False)









