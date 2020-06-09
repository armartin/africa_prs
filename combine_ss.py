import hail as hl
import argparse
from ukbb_pan_ancestry import *


def import_key(ss_filename, ss_keys, ss_name):
    keys = ss_keys.split(',')
    ss = hl.import_table(ss_filename, impute=True, delimiter='\s+',
                         types={keys[1]: hl.tfloat, keys[0]: hl.tstr}, min_partitions=100)
    ss = ss.annotate(**{keys[1]: hl.int(ss[keys[1]])})
    chroms = set(map(str, range(1, 23)))
    ss = ss.filter(hl.literal(chroms).contains(ss[keys[0]]))
    ss = ss.select(locus=hl.locus(hl.str(ss[keys[0]]), ss[keys[1]]), alleles=[ss[keys[2]], ss[keys[3]]],
                     **{'p_' + ss_name: ss[keys[4]], 'beta_' + ss_name: ss[keys[5]]})
    ss = ss.key_by(ss.locus, ss.alleles)
    return ss


def main(args):
    ss = args.ss.split(',')
    chr_pos_ref_alt_p_beta = args.chr_pos_ref_alt_p_beta.split(';')
    ss_names = args.ss_names.split(',')
    sumstats = []

    #  read in each set of sumstats
    for sumstat in range(len(ss)):
        ss_data = import_key(ss[sumstat], chr_pos_ref_alt_p_beta[sumstat], ss_names[sumstat])
        sumstats.append(ss_data)

    ss_joined = sumstats[0]
    for sumstat in range(1, len(ss)):
        ss_joined = ss_joined.join(sumstats[sumstat], 'outer')

    ss_joined = annotate_nearest_gene(ss_joined, add_only_gene_symbols_as_str=True)
    ss_joined = ss_joined.key_by()
    ss_joined = ss_joined.select(chrom=ss_joined.locus.contig, pos=ss_joined.locus.position,
                                 ref=ss_joined.alleles[0], alt=ss_joined.alleles[1],
                                 nearest_genes=ss_joined.nearest_genes,
                                 **ss_joined.row.drop('locus', 'alleles', 'nearest_genes'))

    p_colnames = [x for x in ss_joined.row if x.startswith('p_')]
    ss_filt = ss_joined.filter(hl.sum([hl.is_defined(ss_joined[x]) for x in p_colnames]) > 1)

    ss_filt.export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss', help='comma-separated list of sumstats')
    parser.add_argument('--chr_pos_ref_alt_p_beta', help='separate within dataset by comma, across dataset by ;')
    parser.add_argument('--ss_names', help='comma-separated list of cohort names')
    parser.add_argument('--out')
    args = parser.parse_args()

    main(args)
