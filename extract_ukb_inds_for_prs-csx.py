import hail as hl

target_samples = hl.import_table('gs://apcdr/ukb_holdout/ukb31063.gwas_samples.holdout_and_target.txt', key='s')

contig = 'autosomes'
contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'

ht_variants = hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht')

mt = hl.import_bgen(
        path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
        sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
        entry_fields=['dosage'],
        variants=ht_variants)

mt_target = mt.filter_cols(hl.is_defined(target_samples[mt.s]))  # target

hl.export_plink(mt_target, 'gs://apcdr/ukb_holdout/ukb31063.holdout.target_individuals', ind_id= mt_target.s,
                varid=mt_target.rsid)