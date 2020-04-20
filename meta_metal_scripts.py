"""
take combinations of cohorts and write metal scripts to perform inverse variance-weighted meta-analysis
"""

import argparse
import gzip

def write_cohort_info(out, path, ukb):
    out.write('MARKER\trsid\n')
    out.write('EFFECT\tbeta\n')
    out.write('STDERR\tse\n')
    if ukb:
        #out.write('WEIGHT\tn_complete_samples\n')  # check this
        out.write('ALLELE\tref alt\n')
        out.write('FREQ\tmaf\n')
        out.write('PVAL\tpval\n')
    else:
        #out.write('WEIGHT\tn_complete_samples\n')  # check this
        out.write('ALLELE\tallele0 allele1\n')
        out.write('FREQ\taf\n')
        out.write('PVAL\tp_lrt\n')  # not score


    out.write('PROCESS ' + path + '\n\n')


def main(args):
    ss_inputs = open(args.ss_inputs)
    ss_inputs.readline()
    for line in ss_inputs:
        line = line.strip().split()
        pheno = line[0]
        apcdr = line[1]
        ukb = line[2]
        out = open(args.dirname + '/' + pheno + '.APCDR_UKB.metal.txt', 'w')

        out.write('SCHEME STDERR\n\n')
        write_cohort_info(out, apcdr, False)
        write_cohort_info(out, ukb, True)
        out.write('OUTFILE ' + args.dirname  + '/' + pheno + '_APCDR_UKB.meta .tbl\n')
        out.write('ANALYZE HETEROGENEITY\n')
        out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss_inputs', help='gives path and filename to metal setup file (pheno, sumstat1, sumstat2')
    parser.add_argument('--dirname', help='directory name where output scripts will be written')
    args = parser.parse_args()
    main(args)

