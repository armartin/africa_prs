qimport hail as hl
import argparse
import pandas as pd


def irnt_funct(he: hl.expr.Expression, output_loc: str = 'irnt'):
    ht = he._indices.source
    n_rows = ht.count()
    ht = ht.order_by(he).add_index()
    return ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})


def get_overlapping_phenos(apcdr_ukb, gwas_phenos, gwas_biomarkers, pheno_table, overwrite):
    # get which phenotypes exist in apcdr data
    pheno_gwas = hl.import_table(apcdr_ukb)
    pheno_gwas = {row['pheno_code']: row['ukb_code'] for row in pheno_gwas.collect()}
    incl_whr = False
    if 'WHR' in pheno_gwas:
        del pheno_gwas['WHR']
        incl_whr = True

    if 'EA' in pheno_gwas:
        educ = hl.import_table(
            'gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_final_output/January_2020_plus_pharma_and_updated_codings/phesant_output_multi_ancestry_combined_both_sexes_no_sex_specific_educ_att_years.tsv',
            missing='', impute=True, min_partitions=100, key='userId', types={'userId': hl.tstr, 'EDUC_ATT_CAT_ORD': hl.tint})

    # read ukb data
    ht_phenotypes = hl.import_table(gwas_phenos, force_bgz=True, missing='', impute=True,
                                    min_partitions=100, types={'s': hl.tstr}, key='s')

    phenotype_cols = set(ht_phenotypes.row)
    irnt = []
    raw_phenos = []
    biomarker = []
    for pheno in pheno_gwas.values():
        if pheno + '_irnt' in phenotype_cols:
            irnt.append(pheno)
        elif pheno in phenotype_cols:
            raw_phenos.append(pheno)
        elif pheno + '_0' in phenotype_cols:
            raw_phenos.append(pheno + '_0')
        else:
            biomarker.append(pheno)

    print(pheno_gwas.values())
    if incl_whr:
        irnt.append('whr')

        ht_phenotypes = ht_phenotypes.annotate(whr=ht_phenotypes['48_raw'] / ht_phenotypes['49_raw'])
        ht_phenotypes = irnt_funct(ht_phenotypes.whr, 'whr_irnt').key_by('s')

    if 'EA' in pheno_gwas:
        ht_phenotypes = ht_phenotypes.annotate(ea=educ[ht_phenotypes.key].EDUC_ATT_CAT_ORD)
        raw_phenos.append('ea')
        print('adding EA')

    # now select phenotypes that are in apcdr data
    ht_phenos = ht_phenotypes.select(*[x + '_irnt' for x in irnt] + raw_phenos)
    ht_phenos.show()

    # filter biomarkers to codes
    biomarkers = ['cholesterol_irnt', 'hdl_cholesterol_irnt', 'ldl_irnt', 'triglycerides_irnt',
                  'albumin_irnt', 'alkaline_phosphatase_irnt', 'alanine_aminotransferase_irnt',
                  'aspartate_aminotransferase_irnt', 'direct_bilirubin_irnt',
                  'gamma_glutamyltransferase_irnt', 'glycated_haemoglobin_irnt']
    if gwas_biomarkers:
        ht_biomarkers = hl.read_table(gwas_biomarkers).select(*biomarkers)

        # now join biomarkers with phenotypes
        ht_all_phenos = ht_phenos.join(ht_biomarkers)
        ht_all_phenos.write(pheno_table, overwrite=args.overwrite)
    else:
        ht_phenos.write(pheno_table, overwrite=args.overwrite)


def run_grouped_regressions(mt, ss_output, pheno, pheno_name):
    ht = hl.linear_regression_rows(
        y=[[mt['phenotypes'][y]] for y in pheno],
        x=mt.dosage,
        covariates=[1, *[mt['covariates'][x] for x in list(mt['covariates'].keys())]],
        pass_through=['varid', 'rsid'])

    ht = ht.annotate_globals(phenotypes=pheno) # check this

    ht.write(ss_output + pheno_name + '.ht', overwrite=args.overwrite)


def main(args):
    # get phenotypes that overlap with APCDR dataset
    if args.write_phenos:
        get_overlapping_phenos(args.pheno_ukb_codes, args.gwas_phenos, args.gwas_biomarkers, args.pheno_table,
                               args.overwrite)

    ht_phenos = hl.read_table(args.pheno_table)

    ht_covariates = hl.read_table(args.gwas_covariates)
    ht_variants = hl.read_table(args.gwas_variants)
    ht_samples = hl.import_table(args.gwas_samples, types={'s': hl.tstr}, key='s')

    contig = 'autosomes'
    contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'

    # import ukb bgen files
    mt = hl.import_bgen(
        path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
        sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
        entry_fields=['dosage'],
        variants=ht_variants)

    # add phenotype and covariate info
    mt = mt.annotate_cols(
        phenotypes=ht_phenos[mt.s],
        covariates=ht_covariates[mt.s])

    # filter keeping samples in gwas, i.e. exclude holdout samples
    # mt = mt.filter_cols(ht_samples[mt.s].in_gwas == 'TRUE')
    disc_or_target = 'target'
    if disc_or_target == 'target':
        mt = mt.filter_cols(ht_samples[mt.s].in_gwas == 'FALSE')  # target
    else:
        mt = mt.filter_cols(ht_samples[mt.s].in_gwas == 'TRUE')  # discovery


    # mt.write(args.mt, overwrite=args.overwrite)
    #
    # mt = hl.read_matrix_table(args.mt)

    phenotypes = list(mt['phenotypes'].keys())

    pheno1 = phenotypes[0:10]
    pheno2 = phenotypes[10:20]
    pheno3 = phenotypes[20:30]
    pheno4 = phenotypes[30:len(phenotypes)]
    pheno_leftover = ['whr_irnt', 'glycated_haemoglobin_irnt']

    # run_grouped_regressions(mt, args.holdout_ss_output, pheno1, 'pheno1')
    # run_grouped_regressions(mt, args.holdout_ss_output, pheno2, 'pheno2')
    # run_grouped_regressions(mt, args.holdout_ss_output, pheno3, 'pheno3')
    # run_grouped_regressions(mt, args.holdout_ss_output, pheno4, 'pheno4')
    # run_grouped_regressions(mt, args.holdout_ss_output, pheno_leftover, 'pheno_leftover')
    if disc_or_target == 'target':
        run_grouped_regressions(mt, args.holdout_ss_output, phenotypes, 'target_holdout2')
    else:
        run_grouped_regressions(mt, args.holdout_ss_output, phenotypes, 'discovery')


if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--write_phenos', action='store_true')
    parser.add_argument('--write_mt', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--pheno_table', default='gs://apcdr/ukb_holdout/apcdr_pheno.ht')
    parser.add_argument('--mt', default='gs://apcdr/ukb_holdout/ukb_apcdr_pheno.mt')
    parser.add_argument('--pheno_ukb_codes', default='gs://apcdr/pheno_code_ukb_code.txt')
    parser.add_argument('--gwas_phenos', default='gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz')
    parser.add_argument('--gwas_biomarkers', default='gs://ukb31063/ukb31063.biomarkers_gwas.both_sexes.ht')
    parser.add_argument('--gwas_genos', default='gs://ukb31063/ukb31063.genotype.mt')
    parser.add_argument('--gwas_covariates', default='gs://ukb31063/ukb31063.neale_gwas_covariates.both_sexes.ht')
    parser.add_argument('--gwas_variants', default='gs://ukb31063/ukb31063.neale_gwas_variants.ht')
    parser.add_argument('--gwas_samples', default='gs://apcdr/ukb_holdout/ukb31063.gwas_samples.gwas_vs_holdout.txt')
    parser.add_argument('--holdout_ss_output', default='gs://apcdr/ukb_holdout/ukb31063.gwas_holdout_sumstats_')
    args = parser.parse_args()
    main(args)

