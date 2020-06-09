import hail as hl
import argparse
import numpy as np
# !/opt/conda/default/bin/pip install matplotlib
import matplotlib.pyplot as plt


def import_key(ss_filename, ss_keys):
    ss = hl.import_table(ss_filename, impute=True, delimiter='\s+')
    keys = ss_keys.split(',')
    p = keys[-1]
    ss = ss.annotate(locus=hl.locus(hl.str(ss[keys[0]]), ss[keys[1]]), alleles=[ss[keys[2]], ss[keys[3]]])
    ss = ss.key_by(ss.locus, ss.alleles)
    return ss, p


def main(args):
    ss1, p1 = import_key(args.ss1, args.ss1_chr_pos_ref_alt_p)
    ss2, p2 = import_key(args.ss2, args.ss2_chr_pos_ref_alt_p)
    ss1 = ss1.annotate(ss2=ss2[ss1.key])
    x = (-hl.log10(ss1[p1])).collect()
    y = (-hl.log10(ss1.ss2[p2])).collect()

    fig, ax = plt.subplots()
    plt.xlabel(args.ss1_name)
    plt.ylabel(args.ss2_name)
    plt.title(args.trait)
    ax.scatter(x, y)
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    out_base = args.out.split('/')[-1]
    fig.savefig('/tmp/' + out_base)
    hl.hadoop_copy('file:///tmp/' + out_base, args.out)

    # hl.plot.scatter(-hl.log10(ss1[p1]), -hl.log10(ss1.ss2[p2]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss1', default='gs://apcdr/prs_sumstats_clumps/bbj_page_ukb_10k_eur_holdout_meta/MCHC.meta.bgz')
    parser.add_argument('--ss2')
    parser.add_argument('--ss1_chr_pos_ref_alt_p', default='CHR,BP,A1,A2,P')
    parser.add_argument('--ss2_chr_pos_ref_alt_p', default='chrom,pos,ref,alt,pval')
    parser.add_argument('--ss1_name', default='BBJ+PAGE+UKB')
    parser.add_argument('--ss2_name', default='UKB')
    parser.add_argument('--trait')
    parser.add_argument('--out')
    args = parser.parse_args()

    main(args)
