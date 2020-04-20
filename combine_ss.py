import hail as hl
#from ukb_common import *

pheno_gwas = hl.import_table(f'gs://apcdr/pheno_code_ukb_code.txt')
phenos = pheno_gwas.pheno_code.collect()

alp = hl.import_table('gs://apcdr/apcdr_ukb_10k_eur_holdout_meta/ALP.meta.bgz', impute=True, delimiter='\s+')

row_keys = ['locus', 'alleles']
col_keys = ['pheno', 'pheno_description']
all_hts = list(map(lambda x: hl.read_table(x), all_ht_paths))
rekeyed_hts = pull_out_col_keys(all_hts, row_keys, col_keys)
mt = mwzj_hts_by_tree(rekeyed_hts, temp_dir, col_keys, repartition_final=n_partitions)