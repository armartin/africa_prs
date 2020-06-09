import argparse
import itertools
import gzip
import copy


def make_all_combos(list_to_combine):
    combos = list(itertools.combinations(list_to_combine, 2))
    if len(list_to_combine) > 2:
        for i in range(3, len(list_to_combine) + 1):
            combos.extend(list(itertools.combinations(list_to_combine, i)))
    return combos


def intersect_clumps(clump_dict, window_size, clumps, pop_order):
    header = []
    values = []

    combos = make_all_combos(clump_dict.keys())

    if pop_order:
        combos = make_all_combos(list(zip(clumps, pop_order)))

    for combo in range(len(combos)):

        if pop_order:
            combo_names = [x[1] for x in combos[combo]]
            combo_filenames = [x[0] for x in combos[combo]]
            #print(combo_names)
            header.append('_'.join(combo_names))
            #print(combo_names[combo])
        else:
            combo_filenames = combos[combo]
        clump_dict_subset = {}
        for subset in combo_filenames:
            clump_dict_subset[subset] = clump_dict[subset]

        union_chrs = set()
        for clump in clump_dict_subset.keys():
            union_chrs |= set(clump_dict_subset[clump].keys())

        #  get intersecting chrs
        #  within intersecting chrs, check if each locus is within window_size. start with first in tuple

        intersect_chrs = set()
        clump_filenames = list(clump_dict_subset.keys())
        intersect_chrs |= set(clump_dict_subset[clump_filenames[0]])
        for clump in clump_dict_subset.keys():
            intersect_chrs &= set(clump_dict_subset[clump].keys())

        total_intersect_loci = {}

        for chr in intersect_chrs:
            intersect_snps = []
            intersect_snps.extend(clump_dict_subset[clump_filenames[0]][chr])
            for clump in clump_filenames[1:]:
                current_snps = clump_dict_subset[clump][chr]
                for snp in copy.deepcopy(intersect_snps):
                    # check each of current_snps. if non are within 250 kb, delete snp.
                    keep = False
                    for snp2 in current_snps:
                        if snp2 > snp - window_size and snp2 < snp + window_size:
                            keep = True
                            break
                    if not keep:
                        intersect_snps.remove(snp)
            total_intersect_loci[chr] = intersect_snps
        #print('# total intersect loci: ' + str(sum([len(x) for x in total_intersect_loci.values()])))
        values.append(sum([len(x) for x in total_intersect_loci.values()]))
    return header, values


def independent_clumps(clump_dict, window_size):
    unique_clumps = copy.deepcopy(clump_dict)
    for chr in clump_dict.keys():
        clump_dict[chr].sort()
        last_clump = clump_dict[chr][0]
        if len(clump_dict[chr]) > 1:
            for current_clump in clump_dict[chr][1:len(clump_dict[chr])]:
                if current_clump < last_clump + window_size:
                    unique_clumps[chr].remove(current_clump)
                else:
                    last_clump = current_clump
    # print('# unique loci: ' + str(sum([len(x) for x in unique_clumps.values()])))
    return(unique_clumps)


def make_clump_dict(clump, gwsig):
    """
    Reads a clumped file from plink into a dict
    :param clump:
    :param gwsig:
    :return:
    """
    if clump.endswith('gz'):
        clump_file = gzip.open(clump)
    else:
        clump_file = open(clump)
    header = {k: v for v, k in enumerate(clump_file.readline().strip().split())}
    #clump_dict = {'TOTAL': 0}
    clump_dict = {}
    for line in clump_file:
        line = line.strip().split()
        if float(line[header['P']]) < gwsig:
            if line[header['CHR']] in clump_dict:
                clump_dict[line[header['CHR']]].append(int(line[header['BP']]))
            else:
                clump_dict[line[header['CHR']]] = [int(line[header['BP']])]
            #clump_dict['TOTAL'] += 1
        else:
            break
    return(clump_dict)


def main(args):
    out = open(args.out, 'w')
    header = ['pheno']
    values = [args.pheno]
    clumps = args.clumps.split(',')
    if args.pop_order:
        pop_order = args.pop_order.split(',')
        header.extend(pop_order)
    else:
        header.extend(clumps)
    clump_dicts = {}
    for clump in clumps:
        clump_dict = make_clump_dict(clump, args.gwsig)

        clump_dict_indep = independent_clumps(clump_dict, args.window_size)
        values.append(sum([len(x) for x in clump_dict_indep.values()]))
        clump_dicts[clump] = clump_dict_indep

    intersect_header, intersect_values = intersect_clumps(clump_dicts, args.window_size, clumps, pop_order)
    header.extend(intersect_header)
    values.extend(intersect_values)

    out.write('\t'.join(header) + '\n')
    out.write('\t'.join(map(str, values)) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--clumps', help='comma-separated list of clumped files from different sumstats')
    parser.add_argument('--pop_order', help='list of populations in same order as pops '
                                            '(if not given, alternative just uses clumps)')
    parser.add_argument('--gwsig', default=5e-8, type=float)
    parser.add_argument('--window_size', default=250000, type=int)
    parser.add_argument('--pheno')
    parser.add_argument('--out')
    args = parser.parse_args()

    main(args)
