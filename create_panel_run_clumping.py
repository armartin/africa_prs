import argparse
import subprocess
import random

#check_call
#check_output


def random_pops(analysis, ss_line, ss_info_header_dict, tgp_dict, pop, tmp_dir, append=False):
    """
    Randomly samples the correct number of individuals given a population from 1000 Genomes data
    :param analysis: folder name for each analysis (e.g. meta-analysis) containing sumstats
    :param ss_line: sumstats line info with per-population counts to sample
    :param ss_info_header_dict: dict of header from sumstats to their indices
    :param tgp_dict: 1000 Genomes dict of dicts with super population -> individual ID -> population
    :param pop: population string from EUR, EAS, or AFR
    :param tmp_dir: temporary directory with text file consisting of individuals to keep
    :return:
    """
    num_pop = int(ss_line[ss_info_header_dict['n_' + pop]])
    sampled_inds = random.sample(list(tgp_dict[pop].keys()), num_pop)
    out_file = tmp_dir + analysis + '_' + ss_line[0] + '.keep'
    if not append:
        out = open(out_file, 'w')
    else:
        out = open(out_file, 'a')
    for ind in sampled_inds:
        out.write(tgp_dict[pop][ind] + '\t' + ind + '\n')
    out.close()


def main(args):
    tgp_info = open(args.tgp_info)
    ss_info = open(args.ss_info)

    #  load relevent 1kG reference panel individual info
    tgp_dict = {'EUR': {}, 'EAS': {}, 'AFR': {}, 'AMR': {}}
    for line in tgp_info:
        line = line.strip().split()
        if line[2] in tgp_dict.keys() and line[1]:
            tgp_dict[line[2]][line[0]] = line[1]

    #  load sumstats info
    ss_info_header = ss_info.readline().strip().split()
    ss_info_header_dict = dict(zip(ss_info_header, range(len(ss_info_header))))
    ss_dict = {}
    for line in ss_info:
        line = line.strip().split('\t')
        ss_dict[line[0]] = line

    #  choose individuals for LD reference panels corresponding to each phenotype and meta-analysis
    analysis_folders = ['apcdr_bbj_ukb_10k_eur_holdout_meta', 'apcdr_ukb_10k_eur_holdout_meta',
                        'bbj_ukb_10k_eur_holdout_meta', 'bbj_page_ukb_10k_eur_holdout_meta']
    for pheno in ss_dict:
        #for analysis in analysis_folders:
        for analysis in ['bbj_page_ukb_10k_eur_holdout_meta']:
            if args.pops_to_clump == 'apcdr_bbj_ukb_10k_eur_holdout_meta':
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EAS', args.tmp_dir, append=True)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'AFR', args.tmp_dir, append=True)
            elif args.pops_to_clump == 'apcdr_ukb_10k_eur_holdout_meta':
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'AFR', args.tmp_dir,
                            append=True)
            elif args.pops_to_clump == 'bbj_ukb_10k_eur_holdout_meta':
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EAS', args.tmp_dir,
                            append=True)
            elif args.pops_to_clump == 'bbj_page_ukb_10k_eur_holdout_meta':
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EAS', args.tmp_dir, append=True)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'AFR', args.tmp_dir, append=True)
                random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'AMR', args.tmp_dir, append=True)
            # elif args.pops_to_clump == 'ukb_10k_eur_holdout_meta':
            #     random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)
            # #  this case a little trickier:
            # #  - add folder to sumstats dir
            # #  - use all EAS to clump rather than ss_info
            # elif args.pops_to_clump == 'bbj':
            #     random_pops(analysis, ss_dict[pheno], ss_info_header_dict, tgp_dict, 'EUR', args.tmp_dir)

    #  use plink to make phenotype- and analysis-specific LD reference panels
    #for analysis in analysis_folders:
    for analysis in ['bbj_page_ukb_10k_eur_holdout_meta']:
        print(analysis)
        subprocess.check_call(
            'qsub /humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/ld_panel_clump.sh ' + analysis,
            shell=True
        )


if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--tgp_info', default='/humgen/atgu1/fs03/shared_resources/1kG/integrated/20130502/integrated_call_samples_v3.20130502.ALL.panel')
    parser.add_argument('--ss_info', default='/humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/meta_analyze_ukb_bbj_ugr.txt')
    parser.add_argument('--pops_to_clump', default='apcdr_bbj_ukb_10k_eur_holdout_meta')
    parser.add_argument('--tgp_geno', default='/humgen/atgu1/fs03/shared_resources/1kG/integrated/20130502/ALL.1KG_phase3.20130502.genotypes.maf005')
    parser.add_argument('--tmp_dir', default='/humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/clump_panel/')
    args = parser.parse_args()
    main(args)
