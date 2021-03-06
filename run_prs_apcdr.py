import hail as hl
import pickle
import time
from datetime import datetime
import argparse
from pprint import pprint

def specific_clumps(filename):
    clump = hl.import_table(filename, delimiter='\s+', min_partitions=10,
                            types={'P': hl.tfloat, 'BETA': hl.tfloat})
    clump_dict = clump.aggregate(hl.dict(hl.agg.collect(
        (hl.locus(hl.str(clump.CHR), hl.int(clump.BP)),
        True)
    )), _localize=False)
    return clump_dict


def main(args):
    ########################################################################
    ### initialize
    print('Getting started: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # 1. Read in summary stats data
    # 2. Annotate matrix table with effect sizes for each phenotype
    # 3. Compute PRS for each
    start = time.time()

    pheno_gwas = hl.import_table(f'gs://apcdr/pheno_code_ukb_code.txt')
    pheno_ss = dict([(x.pheno_code, x.ukb_code) for x in pheno_gwas.collect()])
    #pheno_ss = dict([(x.ss_code, x.pheno_code) for x in pheno_gwas.collect()])

    # mt = hl.read_matrix_table('gs://apcdr/prs_sumstats_clumps/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37_subset.mt')
    mt = hl.read_matrix_table('gs://apcdr/dosage_bgen/apcdr.mt')
    ss_keys = dict(zip(['CHR', 'POS', 'REF', 'ALT', 'P', 'BETA'], args.chr_pos_ref_alt_p_beta.split(',')))

    for pheno in list(pheno_ss.keys()):
    #for pheno in ['WHR']:
        print('Pheno: ' + pheno + ', Time: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        suffix_replace = args.ss_suffix.split('.')
        suffix_replace[-2] = 'clumped'
        suffix_replace = '.'.join(suffix_replace)
        if hl.hadoop_exists(args.ss_clump_prefix + pheno + suffix_replace):
            ss_path = args.ss_clump_prefix + pheno + args.ss_suffix
            clump_path = args.ss_clump_prefix + pheno + suffix_replace
        elif hl.hadoop_exists(args.ss_clump_prefix + pheno_ss[pheno] + suffix_replace):
            ss_path = args.ss_clump_prefix + pheno_ss[pheno] + args.ss_suffix
            clump_path = args.ss_clump_prefix + pheno_ss[pheno] + suffix_replace
        else:
            continue

        ss = hl.import_table(ss_path, impute=True, delimiter='\s+', min_partitions=1000)
        ss = ss.annotate(locus = hl.locus(hl.str(ss[ss_keys['CHR']]), ss[ss_keys['POS']]), alleles=[ss[ss_keys['REF']], ss[ss_keys['ALT']]])
        ss = ss.key_by(ss.locus, ss.alleles)

        ## Read in summary statistics and true phenotypes
        mt_annot = mt.annotate_rows(ss=ss[mt.locus, mt.alleles]) # come back to this
        # ht_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.gwas_vs_holdout.txt',
        #                              types={'s': hl.tstr}, key='s')
        # ht_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.holdout_and_target.txt',
        #                              types={'s': hl.tstr}, key='s')
        # # mt_annot = mt_annot.filter_cols(hl.or_else(ht_samples[mt_annot.s].in_gwas != 'TRUE', True))
        # mt_annot = mt_annot.filter_cols(hl.is_defined(ht_samples[mt_annot.s]))
        # # print(mt.count()) # 13364303, 136265)

        print('Starting ' + pheno + ': ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5,
                 's10': 1}

        pheno_clump = specific_clumps(clump_path)

        mt_annot = mt_annot.filter_rows(pheno_clump.get(mt_annot.locus, False))
        # print(mt.count())

        annot_expr = {
            k: hl.agg.sum(hl.float(mt_annot.ss[ss_keys['BETA']]) * mt_annot.dosage * hl.int(
                mt_annot.ss[ss_keys['P']] < v))
            for k, v in p_max.items()}

        mt_annot = mt_annot.annotate_cols(**annot_expr)

        ht_out = mt_annot.cols()
        #ht_out.describe()
        #covs = hl.read_table('gs://apcdr/ukb_holdout/uk_round2_allSamples_phenos_phesant.ht').select('age', 'sex')  # added
        # need to add in PCs
        #ht_out = ht_out.annotate(**covs[ht_out.key])
        ht_comb = ht_out.select(*p_max.keys(),
                               age=ht_out.phenotypes.age,
                               sex=ht_out.phenotypes.sex,
                               pheno=ht_out.phenotypes[pheno])

        output_location = args.ss_clump_prefix + pheno + '_apcdr_PRS'
        #ht_comb.describe()
        #ht_comb.write(output_location + '.ht', overwrite=args.overwrite)
        #ht_comb = hl.read_table(output_location + '.ht')
        ht_comb.export(output_location + '.txt.bgz')

    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss_clump_prefix', default='gs://apcdr/prs_sumstats_clump/apcdr_ukb_10k_eur_holdout_meta/')
    parser.add_argument('--ss_suffix', default='.meta.bgz')
    parser.add_argument('--chr_pos_ref_alt_p_beta', default='CHR,POS,A1,A2,P,BETA')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    hl.init(log='/prs.log')
    main(args)