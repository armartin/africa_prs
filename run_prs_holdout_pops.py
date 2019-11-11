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
    )))
    return hl.literal(clump_dict)


def main(args):
    ########################################################################
    ### initialize
    print('Getting started: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    apcdr_phenos = {'21001': '21001_irnt', '21002': '21002_irnt', '30000': '30000_irnt', '30010': '30010_irnt',
                    '30020': '30020_irnt', '30030': '30030_irnt', '30040': '30040_irnt', '30050': '30050_irnt',
                    '30060': '30060_irnt', '30070': '30070_irnt', '30080': '30080_irnt', '30100': '30100_irnt',
                    '30120': '30120_irnt', '30130': '30130_irnt', '30140': '30140_irnt', '30150': '30150',
                    '30160': '30160', '30180': '30180_irnt', '30190': '30190_irnt', '30200': '30200_irnt',
                    '30210': '30210_irnt', '30220': '30220_irnt', '4079': '4079_irnt', '4080': '4080_irnt',
                    '48': '48_irnt', '49': '49_irnt', '50': '50_irnt', '30690': 'cholesterol_irnt',
                    '30760': 'hdl_cholesterol_irnt', '30780': 'ldl_irnt', '30870': 'triglycerides_irnt',
                    '30600': 'albumin_irnt', '30610': 'alkaline_phosphatase_irnt',
                    '30620': 'alanine_aminotransferase_irnt', '30650': 'aspartate_aminotransferase_irnt',
                    '30840': 'direct_bilirubin_irnt', '30730': 'gamma_glutamyltransferase_irnt'}

    start = time.time()
    ss = hl.read_table('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37.ht')  # added

    if args.mt_comb:
        contig = 'autosomes'
        contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'

        start = time.time()
        # large block size because we read very little data (due to filtering & ignoring genotypes)
        hl.init(branching_factor=10, min_block_size=2000)

        mt_all = hl.import_bgen(
            path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
            sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
            entry_fields=['dosage'],
            variants=hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht'))

        ## Read in summary statistics and true phenotypes
        all_phenos = hl.read_table('gs://apcdr/ukb_holdout/uk_round2_allSamples_phenos_phesant.ht') # added

        mt_all = mt_all.annotate_cols(phenotypes = all_phenos[mt_all.s])
        mt_all = mt_all.annotate_rows(ss=ss[mt_all.locus, mt_all.alleles])
        mt_all.write('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37.mt', overwrite=args.overwrite)

    mt_all = hl.read_matrix_table('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_pheno37.mt')

    phenotypes = ss.phenotypes.take(50)[0]

    print('Read in bgen, sumstat, pheno info: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for pheno in apcdr_phenos.keys():
        print('Starting ' + pheno + ': ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1}

        pheno_clump = specific_clumps('gs://apcdr/ukb_holdout/ukb31063.gwas_holdout.' + pheno + '_v3.clumped')

        mt = mt_all.filter_rows(pheno_clump.get(mt_all.locus, False))
        # print(mt.count())

        annot_expr = {
            k: hl.agg.sum(mt.ss.beta[phenotypes.index(apcdr_phenos[pheno])][0] * mt.dosage * hl.int(mt.ss.p_value[phenotypes.index(apcdr_phenos[pheno])][0] < v))
            for k, v in p_max.items()}

        mt = mt.annotate_cols(**annot_expr)

        ht_out = mt.cols()
        #ht_out.describe()
        #ht_out['phenotypes'].show()
        ht_prs = ht_out.select(*p_max.keys())
        ht_covars = ht_out.select(age=ht_out.phenotypes.age,
                                  sex=ht_out.phenotypes.sex,
                                  pheno=ht_out.phenotypes[pheno])
        #ht_covars.show()
        ht_comb = ht_prs.join(ht_covars)

        output_location = 'gs://apcdr/ukb_holdout/ukb31063.gwas_holdout.' + pheno + '_PRS'
        print(ht_comb.describe())
        ht_comb.write(output_location + '.ht')
        ht_comb = hl.read_table(output_location + '.ht')
        ht_comb.export(output_location + '.txt.bgz')

    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--mt_comb', action='store_true')
    args = parser.parse_args()

    hl.init(log='/prs.log')
    main(args)