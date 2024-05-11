import sys
import pandas as pd

header = """##fileformat=VCFv4.2
##fileDate=20231016
##reference=GRCh38/hg38
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""


if __name__ == '__main__':
    vars_df = pd.read_csv(sys.argv[1], sep='\t')
    with open(sys.argv[2], 'w') as f: #vcf
        f.write(header)
        for i,r in vars_df.iterrows():
            f.write('\t'.join(map(str, [r['variant'].split('_')[0].strip('chr'), r['variant'].split('_')[1], '.', r['variant'].split('_')[-2], r['variant'].split('_')[-1], '.', '.', '.'])))
            f.write('\n')