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
    pheno_ss = dict([(x.ukb_code, x.pheno_code) for x in pheno_gwas.collect()])
    #pheno_ss = dict([(x.ss_code, x.pheno_code) for x in pheno_gwas.collect()])

    mt = hl.read_matrix_table('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37_subset.mt')


    for pheno in list(pheno_ss.values()):
    #for pheno in ['WHR']:
        print('Pheno: ' + pheno + ', Time: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        if hl.hadoop_exists(args.ss_clump_prefix + pheno + '.meta.bgz'):
            ss = hl.import_table(args.ss_clump_prefix + pheno + '.meta.bgz', impute=True, delimiter='\s+',
                                 min_partitions=1000)
            if 'A1' in list(ss.row):
                ss = ss.annotate(locus = hl.locus(hl.str(ss.CHR), ss.BP), alleles=[ss.A1, ss.A2])
            else:
                ss = ss.annotate(locus=hl.locus(hl.str(ss.CHR), ss.POS), alleles=[ss.REF, ss.ALT])
            ss = ss.key_by(ss.locus, ss.alleles)

            ## Read in summary statistics and true phenotypes
            #all_phenos = hl.read_table('gs://apcdr/ukb_holdout/apcdr_pheno.ht')
            all_phenos = hl.read_table('gs://apcdr/ukb_holdout/uk_round2_allSamples_phenos_phesant2.ht')
            all_phenos.describe()
            all_phenos = all_phenos.transmute(**{v: all_phenos[k] for k, v in pheno_ss.items()} )

            mt_annot = mt.annotate_cols(phenotypes = all_phenos[mt.s])
            mt_annot = mt_annot.annotate_rows(ss=ss[mt_annot.locus, mt_annot.alleles]) # come back to this
            # ht_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.gwas_vs_holdout.txt',
            #                              types={'s': hl.tstr}, key='s')
            ht_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.holdout_and_target.txt',
                                         types={'s': hl.tstr}, key='s')
            # mt_annot = mt_annot.filter_cols(hl.or_else(ht_samples[mt_annot.s].in_gwas != 'TRUE', True))
            mt_annot = mt_annot.filter_cols(hl.is_defined(ht_samples[mt_annot.s]))
            # print(mt.count()) # 13364303, 136265)

            print('Starting ' + pheno + ': ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

            p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5,
                     's10': 1}

            pheno_clump = specific_clumps(args.ss_clump_prefix + pheno + '.clumped.bgz')

            mt_annot = mt_annot.filter_rows(pheno_clump.get(mt_annot.locus, False))
            # print(mt.count())

            annot_expr = {
                k: hl.agg.sum(hl.float(mt_annot.ss.BETA) * mt_annot.dosage * hl.int(
                    mt_annot.ss.P < v))
                for k, v in p_max.items()}

            mt_annot = mt_annot.annotate_cols(**annot_expr)

            ht_out = mt_annot.cols()
            covs = hl.read_table('gs://apcdr/ukb_holdout/uk_round2_allSamples_phenos_phesant.ht').select('age', 'sex')  # added
            ht_out = ht_out.annotate(**covs[ht_out.key])
            ht_comb = ht_out.select(*p_max.keys(),
                                   age=ht_out.age,
                                   sex=ht_out.sex,
                                   pheno=all_phenos[ht_out.key][pheno])

            output_location = args.ss_clump_prefix + pheno + '_PRS'
            #ht_comb.describe()
            ht_comb.write(output_location + '.ht', overwrite=args.overwrite)
            ht_comb = hl.read_table(output_location + '.ht')
            ht_comb.export(output_location + '.txt.bgz')

    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss_clump_prefix', default='gs://apcdr/prs_sumstats_clump/apcdr_ukb_10k_eur_holdout_meta/')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    hl.init(log='/prs.log')
    main(args)