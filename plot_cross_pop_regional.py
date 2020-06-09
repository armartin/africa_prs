import hail as hl
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
import seaborn as sns


def import_key(ss_filename, ss_keys, clump_name):
    keys = ss_keys.split(',')
    ss = hl.import_table(ss_filename, impute=True, delimiter='\s+',
                         types={keys[1]: hl.tfloat, keys[0]: hl.tstr}, min_partitions=100)
    clump = hl.import_table(clump_name, delimiter='\s+', min_partitions=10,
                            types={'P': hl.tfloat, 'CHR': hl.tstr, 'BP': hl.tint})
    clump = clump.key_by(locus=hl.locus(clump.CHR, clump.BP))
    clump = clump.filter(clump.P < 5e-8)
    ss = ss.annotate(**{keys[1]: hl.int(ss[keys[1]])})
    chroms = set(map(str, range(1, 23)))
    ss = ss.filter(hl.literal(chroms).contains(ss[keys[0]]))
    ss = ss.annotate(locus=hl.locus(hl.str(ss[keys[0]]), ss[keys[1]]), alleles=[ss[keys[2]], ss[keys[3]]])
    ss = ss.key_by(ss.locus)
    ss = ss.annotate(clump=hl.is_defined(clump[ss.key]))
    ss = ss.key_by(ss.locus, ss.alleles)
    p = keys[-1]
    return ss, p


def main(args):
    #  set up tracking variables corresponding to each dataset
    sumstats = []
    p = []
    ss = args.ss.split(',')
    if args.clumps:
        clumps = args.clumps.split(',')
    ss_names = args.ss_names.split(',')
    chr_pos_ref_alt_p = args.chr_pos_ref_alt_p.split(';')

    #  read in each set of sumstats
    for sumstat in range(len(ss)):
        ss_data, p_data = import_key(ss[sumstat], chr_pos_ref_alt_p[sumstat], clumps[sumstat])
        sumstats.append(ss_data)
        p.append(p_data)

    #  join across datasets, filter to target region
    ss_joined = sumstats[0]
    for sumstat in range(1, len(ss)):
        annot_val = 'ss' + str(sumstat)
        ss_joined = ss_joined.annotate(**{annot_val: sumstats[sumstat][ss_joined.key]})
    ss_joined_filt = ss_joined.filter(hl.parse_locus_interval(args.region).contains(ss_joined.locus))
    ss_to_plot = ss_joined_filt.to_pandas()

    # set up figure and plot
    print('Arial')
    sns.set(font='Arial')
    sns.set_style('white')
    sns.despine()
    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(hspace=0.5)
    plt.tight_layout()
    spacing = 0.04
    fig.text(0.5, spacing, 'Position', ha='center')
    fig.text(spacing, 0.5, '-log10(p)', va='center', rotation='vertical')
    fig.text(0.5, 1 - spacing, args.trait, ha='center')
    # tips = sns.load_dataset("tips")

    for sumstat in range(len(ss)):
        ax = fig.add_subplot(len(ss), 1, sumstat+1)#, sharex='col')
        if sumstat > 0:
            current_p = 'ss' + str(sumstat) + '.' + p[sumstat]
            pos = 'ss' + str(sumstat) + '.' + chr_pos_ref_alt_p[sumstat].split(',')[1]
        else:
            current_p = p[sumstat]
            pos = chr_pos_ref_alt_p[sumstat].split(',')[1]
        ss_current_plot = ss_to_plot.filter(items=[pos, current_p, 'clump'], axis=1)
        ss_current_plot.dropna(inplace=True)
        sns.scatterplot(ss_current_plot[pos],
                        -np.log10(ss_current_plot[current_p]),
                        hue=ss_current_plot.clump,
                        linewidth=0, ax=ax)
        plt.xlabel("")
        plt.ylabel("")
        plt.title(ss_names[sumstat])

    # write plot out
    # gs.tight_layout(fig)
    out_base = args.out.split('/')[-1]
    fig.savefig('/tmp/' + out_base)
    hl.hadoop_copy('file:///tmp/' + out_base, args.out)


    # #  subset to target region and set up plot
    # x = (-hl.log10(ss1[p1])).collect()
    # y = (-hl.log10(ss1.ss2[p2])).collect()
    #
    # fig, ax = plt.subplots()
    # plt.xlabel(args.ss1_name)
    # plt.ylabel(args.ss2_name)
    # plt.title(args.trait)
    # ax.scatter(x, y)
    # lims = [
    #     np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    #     np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    # ]
    # ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    #
    # fig.savefig('/tmp/' + out_base)
    # hl.hadoop_copy('file:///tmp/' + out_base, args.out)

    # hl.plot.scatter(-hl.log10(ss1[p1]), -hl.log10(ss1.ss2[p2]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss', help='comma-separated list of sumstats')
    parser.add_argument('--clumps', help='comma-separated list of clumps, same order')
    parser.add_argument('--chr_pos_ref_alt_p', help='separate within dataset by comma, across dataset by ;')
    parser.add_argument('--ss_names', help='comma-separated list of cohort names')
    parser.add_argument('--trait')
    parser.add_argument('--region', help='e.g. 1:100-1000')
    parser.add_argument('--out')
    args = parser.parse_args()

    main(args)
