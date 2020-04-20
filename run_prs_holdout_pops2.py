import hail as hl
import pickle
import time
from datetime import datetime
import argparse
from pprint import pprint

def specific_clumps(filename):
    clump = hl.import_table(filename, delimiter='\s+', min_partitions=10, types={'P': hl.tfloat})
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
    phenos = dict([(x.pheno_code, x.ukb_code) for x in pheno_gwas.collect()])
    contig = 'autosomes'
    contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'

    # hl.init(branching_factor=10, min_block_size=2000)
    # mt = hl.import_bgen(
    #     path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
    #     sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
    #     entry_fields=['dosage'],
    #     variants=hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht'))
    mt = hl.read_matrix_table('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37_subset.mt')

    for pheno in phenos:
        print('Pheno: ' + pheno + ', Time: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        if hl.hadoop_exists('gs://apcdr/apcdr_ukb_10k_eur_holdout_meta/' + pheno + '.meta.bgz'):
            ss = hl.import_table('gs://apcdr/apcdr_ukb_10k_eur_holdout_meta/' + pheno + '.meta.bgz', impute=True, delimiter='\s+')
            ss = ss.annotate(locus = hl.locus(hl.str(ss.CHR), ss.BP), alleles=[ss.A1, ss.A2])
            ss = ss.key_by(ss.locus, ss.alleles)

            ## Read in summary statistics and true phenotypes
            all_phenos = hl.read_table('gs://apcdr/ukb_holdout/uk_round2_allSamples_phenos_phesant.ht') # added

            mt_annot = mt.annotate_cols(phenotypes = all_phenos[mt.s])
            mt_annot = mt_annot.annotate_rows(ss=ss[mt_annot.locus, mt_annot.alleles]) # come back to this
            ht_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.gwas_vs_holdout.txt',
                                         types={'s': hl.tstr}, key='s')
            mt_annot = mt_annot.filter_cols(hl.or_else(ht_samples[mt_annot.s].in_gwas != 'TRUE', True))
            # print(mt.count()) # 13364303, 136265)

            print('Starting ' + pheno + ': ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

            p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5,
                     's10': 1}

            pheno_clump = specific_clumps('gs://apcdr/apcdr_ukb_10k_eur_holdout_meta/' + pheno + '_clump.clumped.bgz')

            mt_annot = mt_annot.filter_rows(pheno_clump.get(mt_annot.locus, False))
            # print(mt.count())

            annot_expr = {
                k: hl.agg.sum(mt_annot.ss.BETA * mt_annot.dosage * hl.int(
                    mt_annot.ss.P < v))
                for k, v in p_max.items()}

            mt_annot = mt_annot.annotate_cols(**annot_expr)

            ht_out = mt_annot.cols()
            # ht_out.describe()
            # ht_out['phenotypes'].show()
            ht_prs = ht_out.select(*p_max.keys())
            ht_covars = ht_out.select(age=ht_out.phenotypes.age,
                                      sex=ht_out.phenotypes.sex,
                                      pheno=ht_out.phenotypes[phenos[pheno]])
            # ht_covars.show()
            ht_comb = ht_prs.join(ht_covars)
            output_location = 'gs://apcdr/apcdr_ukb_10k_eur_holdout_meta/' + pheno + '_PRS'
            print(ht_comb.describe())
            ht_comb.write(output_location + '.ht', overwrite=args.overwrite)
            ht_comb = hl.read_table(output_location + '.ht')
            ht_comb.export(output_location + '.txt.bgz')

    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    hl.init(log='/prs.log')
    main(args)