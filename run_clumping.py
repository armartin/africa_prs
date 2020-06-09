import argparse
import subprocess


def main(args):
    ss_info = open(args.ss_info)
    ss_info_header = {k: v for v, k in enumerate(ss_info.readline().strip().split('\t'))}
    relevant_path_col = ss_info_header[args.path_colname] + 1
    cmd = 'qsub /home/unix/armartin/atgu/armartin/ginger/apcdr/sumstats/run_clumping2.sh ' + \
        ' '.join([args.ss_info, str(relevant_path_col), args.ref_panel, args.pval_colname, args.snp_colname, args.out])
    print(cmd)
    subprocess.check_call(cmd, shell=True)



if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_panel', default='/home/unix/armartin/atgu/shared_resources/1kG/integrated/20130502/EUR.1KG_phase3.20130502.genotypes.maf005')
    parser.add_argument('--ss_info', default='/humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/meta_analyze_ukb_bbj_ugr.txt')
    parser.add_argument('--path_colname', default='apcdr_file')
    parser.add_argument('--pval_colname', default='P')
    parser.add_argument('--snp_colname', default='SNP')
    parser.add_argument('--out')
    args = parser.parse_args()
    main(args)
