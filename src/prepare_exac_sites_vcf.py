#!/usr/bin/env python

__author__ = 'konradjk'
# With special thanks to Daniel Birnbaum (histograms) and Eric Minikel (distance)

import argparse
import gzip
import pipes
import sys
import numpy
import re
from collections import Counter, defaultdict
import scipy.stats

metrics = ['DP', 'GQ']
bins = range(0, 101, 5)
mid_string = '|'.join(map(str, [float(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]))
bins[-1] = 1000

all_bins = {}
all_bins['DP'] = bins
all_bins['GQ'] = bins

all_mids = {}
all_mids['DP'] = mid_string
all_mids['GQ'] = mid_string


print >> sys.stderr, "Pre-calculating depth cutoffs..."
depth_limit = 250
depths = range(depth_limit + 1)
p_ad_depth = defaultdict(dict)
for depth in depths:
    for x in depths:
        if x <= depth:
            p_ad_depth[depth][x] = -scipy.log(scipy.stats.binom_test(x, depth))
print >> sys.stderr, "Done!"


def main(args):
    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)

    pca_data = read_pcs(args.pca)
    sex_data = read_sex(args.sex)
    consanguineous_samples = read_consanguineous_samples(args.consanguineous)

    # Opening output files
    if not args.output.endswith('.gz'): args.output += '.gz'
    pipe = pipes.Template()
    pipe.append('bgzip -c /dev/stdin', '--')
    g = pipe.open(args.output, 'w')

    header = None
    for line in f:
        line = line.strip()

        # Reading and writing header lines
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                for metric in metrics:
                    print >> g, '##INFO=<ID=%s_HIST,Number=R,Type=String,Description="Histogram for %s; Mids: %s">' % (metric, metric, all_mids[metric])
                print >> g, '##INFO=<ID=DOUBLETON_DIST,Number=A,Type=String,Description="Euclidean distance of carriers of doubletons">'

                print >> g, '##INFO=<ID=AC_MALE,Number=A,Type=String,Description="Allele count among males">'
                print >> g, '##INFO=<ID=AC_FEMALE,Number=A,Type=String,Description="Allele count among females">'

                print >> g, '##INFO=<ID=AN_MALE,Number=1,Type=String,Description="Allele number among males">'
                print >> g, '##INFO=<ID=AN_FEMALE,Number=1,Type=String,Description="Allele number among females">'

                print >> g, '##INFO=<ID=AC_CONSANGUINEOUS,Number=A,Type=String,Description="Allele count among individuals with F > 0.05">'
                print >> g, '##INFO=<ID=AN_CONSANGUINEOUS,Number=1,Type=String,Description="Allele number among individuals with F > 0.05">'
                print >> g, '##INFO=<ID=Hom_CONSANGUINEOUS,Number=A,Type=String,Description="Homozygote count among individuals with F > 0.05">'

                header_list = line.split('\t')
                g.write('\t'.join(header_list[:8]) + '\n')
                header_list = [x.replace('#', '').replace(' ', '_') for x in header_list]
                header = dict([(x.replace('#', '').replace(' ', '_'), i) for i, x in enumerate(header_list)])
            else:
                # Edits for VCF header
                if line.startswith('##INFO=<ID=AC_') or line.startswith('##INFO=<ID=Hom_'):
                    line = line.replace('Number=1', 'Number=A').replace('Type=String', 'Type=Integer')
                elif line.startswith('##INFO=<ID=Het_'):
                    line = line.replace('Number=A', 'Number=.')
                elif line == '##fileformat=VCFv4.1':
                    line = '##fileformat=VCFv4.2'
                g.write(line + '\n')
            continue

        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)

        fields = line.split('\t')
        alt_alleles = fields[header['ALT']].split(',')
        alts = len(alt_alleles)
        # Pull out annotation info from INFO and ALT fields
        new_info = fields[header['INFO']].rstrip(';')

        # Pre-computing histograms
        data_list, ad_means, ad_stdevs = get_histograms_for_variant(fields, metrics, all_bins, alts=alts)
        for i, metric in enumerate(metrics):
            hists = []
            for j in range(alts + 1):
                hist = data_list[i*(alts+1)+j]
                hists.append('|'.join(map(str, hist)))
                new_info += ';%s_HIST=%s' % (metric, ','.join(hists))

        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[header['INFO']])])
        acs = info_field['AC_Adj'].split(',')
        homs = info_field['AC_Hom'].split(',')

        if fields[header['FILTER']] == 'PASS':
            if not sum(map(int, info_field['AC_Adj'].split(','))):
                fields[header['FILTER']] = 'AC_Adj0_Filter'
            elif 'InbreedingCoeff' in info_field and float(info_field['InbreedingCoeff']) <= -0.2:
                fields[header['FILTER']] = 'InbreedingCoeff_Filter'

        doubleton_dists = ['.']*alts

        ac_male = ['.']*alts
        ac_female = ['.']*alts
        ac_consang = ['.']*alts
        hom_consang = ['.']*alts

        all_samples = get_sample_info(fields)
        if hemizygous_x(fields):
            all_samples = dict([(sample, gt) for sample, gt in all_samples.items() if sex_data[header_list[sample]] == 'Female' or len(set(gt)) == 1])
        elif fields[header['CHROM']] == 'Y':
            all_samples = dict([(sample, gt) for sample, gt in all_samples.items() if sex_data[header_list[sample]] == 'Male'])
        # print all_samples
        variant_sex_data = Counter([sex_data[header_list[sample]] for sample in all_samples])
        an_male = variant_sex_data['Male']*2 if not hemizygous_segment(fields) else variant_sex_data['Male']
        an_female = variant_sex_data['Female']*2

        sample_names = set([header_list[sample] for sample in all_samples])
        an_consang = len(set(consanguineous_samples.keys()).intersection(sample_names))*2

        for i, alt in enumerate(alt_alleles):
            allele_num = str(i + 1)
            if acs[i] == '0':
                continue
            samples = dict([(sample, gt) for sample, gt in all_samples.items() if allele_num in gt])
            # Calculate doubleton euclidean distance
            if acs[i] == '2' and homs[i] == '0':
                if len(samples) != 2:
                    print >> sys.stderr, 'Variant %s seems to be AC_Adj = 2, but %s samples found with allele' % ('-'.join([fields[header['CHROM']], fields[header['POS']], fields[header['REF']], alt]), len(samples))
                else:
                    doubleton_samples = samples.keys()
                    if header_list[doubleton_samples[0]] in pca_data and header_list[doubleton_samples[1]] in pca_data:
                        doubleton_dists[i] = euclid_dist(pca_data[header_list[doubleton_samples[0]]], pca_data[header_list[doubleton_samples[1]]])

            # Add male and female allele counts
            ac_male[i] = sum([Counter(gt)[allele_num] for sample, gt in samples.items() if sex_data[header_list[sample]] == 'Male'])
            if hemizygous_segment(fields): ac_male[i] /= 2  # Males will be labelled as homozygous (filtered previously) on non-PAR X/Y
            ac_female[i] = sum([Counter(gt)[allele_num] for sample, gt in samples.items() if sex_data[header_list[sample]] == 'Female'])

            # Get consanguineous counts
            ac_consang[i] = sum([Counter(gt)[allele_num] for sample, gt in samples.items() if header_list[sample] in consanguineous_samples])
            hom_consang[i] = sum([Counter(gt)[allele_num] == 2 for sample, gt in samples.items() if header_list[sample] in consanguineous_samples])

        # Write results
        new_info += ';DOUBLETON_DIST=%s' % (','.join(map(str, doubleton_dists)))
        new_info += ';AC_MALE=%s' % (','.join(map(str, ac_male)))
        new_info += ';AC_FEMALE=%s' % (','.join(map(str, ac_female)))
        new_info += ';AN_MALE=%s;AN_FEMALE=%s' % (an_male, an_female)
        new_info += ';AC_CONSANGUINEOUS=%s;AN_CONSANGUINEOUS=%s;HOM_CONSANGUINEOUS=%s' % (','.join(map(str, ac_consang)), an_consang, ','.join(map(str, hom_consang)))
        fields[header['INFO']] = new_info
        g.write('\t'.join(fields[:8]) + '\n')

    f.close()
    g.close()


# Adapted from Daniel Birnbaum's histogram script
def get_histograms_for_variant(fields, metrics, bins, alts=0):
    indices = [None]*len(metrics)*(alts+1)
    distrs = [[] for _ in range(len(metrics)*(alts+1))]
    output = [None]*len(metrics)*(alts+1)
    for i, metric in enumerate(metrics):
        try:
            if metric.startswith('AD'): metric = 'AD'
            ind = fields[8].split(':').index(metric)
            for j in range(alts + 1):
                indices[i*(alts+1)+j] = ind
        except Exception, e:
            pass
    # get distribution for metric
    for sample in fields[9:]:
        # This is only DP/GQ for now
        sample_info = sample.split(':')
        if sample_info[0] == './.': continue
        gts = map(int, sample_info[0].split('/'))
        gt_counts = Counter(gts)
        for i, metric in enumerate(metrics):
            if indices[i*(alts+1)+j] < len(sample_info) and sample_info[indices[i*(alts+1)+j]] != '.':
                for j in range(alts + 1):
                    if not j or j in gts:
                        datum = sample_info[indices[i*(alts+1)+j]]
                        if metric.startswith('AD'):
                            if gt_counts[j] == 1 and ',' in datum:
                                distrs[i*(alts+1)+j].append(map(int, datum.split(',')))
                        else:
                            distrs[i*(alts+1)+j].append(datum)
    ad_means = []
    ad_stdevs = []
    for i, metric in enumerate(metrics):
        for j in range(alts + 1):
            if metric in ['DP', 'GQ']:
                data_to_hist = map(int, distrs[i*(alts+1)+j])
            elif metric.startswith('AD'):
                ad_sum = [(x[j], sum(x)) for x in distrs[i*(alts+1)+j] if sum(x) > 0]
                if metric == 'AD_FRACTION':
                    data_to_hist = [float(x[0])/x[1] for x in ad_sum]
                    ad_means.append(numpy.mean(data_to_hist))
                    ad_stdevs.append(numpy.std(data_to_hist))
                elif metric == 'AD_PVALUE':
                    # Fixing to max(depths)
                    data_to_hist = []
                    for ad, total in ad_sum:
                        if total > depth_limit:
                            ad = int(round(float(ad)*depth_limit/total))
                            total = depth_limit
                        data_to_hist.append(p_ad_depth[total][ad])
            hist, _ = numpy.histogram(data_to_hist, bins=bins[metric])
            output[i*(alts+1)+j] = map(str, hist)
    return output, ad_means, ad_stdevs


# Adapted from Eric Minikel's Diversity Score scripts
def read_pcs(path, n=9):
    """
    Read principal components from a CSV file at the specified path.
    First column is sample id, next n are principal components. Additional
    columns may be present but will be ignored.
    """
    pcs = {}
    myopen = gzip.open if path.endswith('.gz') else open
    with myopen(path) as inf:
        inf.readline().strip().split(',')
        for line in inf:
            cols = line.strip().split(',')
            sampleid = cols[0]
            samplepcs = [float(col) for col in cols[1:(n+1)]]
            pcs[sampleid] = samplepcs
    return pcs


def read_sex(path):
    """
    Read principal components from a CSV file at the specified path.
    First column is sample id, next n are principal components. Additional
    columns may be present but will be ignored.
    """
    sexes = {}
    myopen = gzip.open if path.endswith('.gz') else open
    with myopen(path) as inf:
        for line in inf:
            cols = line.strip().split('\t')
            sexes[cols[0].replace(' ', '_')] = cols[2]
    return sexes


def read_consanguineous_samples(path, cutoff=0.05):
    """
    Read inbreeding coefficients from a TSV file at the specified path.
    Second column is sample id, 6th column is F coefficient. From PLINK:
    FID, IID, O(HOM), E(HOM), N(NM), F
    Additional columns may be present but will be ignored.
    """
    consanguineous_samples = {}
    myopen = gzip.open if path.endswith('.gz') else open
    with myopen(path) as inf:
        _ = inf.readline()
        for line in inf:
            cols = line.strip().split()
            if float(cols[5]) > cutoff:
                consanguineous_samples[cols[1]] = True
    return consanguineous_samples


def euclid_dist(coords1, coords2, weights=None):
    """
    Given two equal-length lists of coordinates in multi-dimensional space,
    return the Euclidean distance between the two points.
    """
    assert len(coords1) == len(coords2), "Coordinate vectors differ in length"
    squared_diffs = [(coords1[i] - coords2[i])**2 for i in range(0,len(coords1))]
    if weights is not None:
        assert len(weights) == len(squared_diffs), "Weight vector is different length than coordinate vectors"
        squared_diffs = [weights[i]*squared_diffs[i] for i in range(0,len(weights))]
    euclidean_distance = sum(squared_diffs)**.5
    return euclidean_distance


def get_sample_info(fields, depth=10, gq=20):
    samples = {}
    format = dict(zip(fields[8].split(':'), range(len(fields[8].split(':')))))
    for i, x in enumerate(fields[9:]):
        gt = x.split(':')
        try:
            gts = gt[0].split('/')
            if '.' not in gts and int(gt[format['DP']]) >= depth and int(gt[format['GQ']]) >= gq:
                samples[i + 9] = gts
        except ValueError, e:
            pass
    return samples


def hemizygous_segment(fields):
    return fields[0] == 'Y' or hemizygous_x(fields)


def hemizygous_x(fields):
    pos = int(fields[1])
    return fields[0] == 'X' and not ((60001 <= pos <= 2699520) or (154931044 <= pos <= 155260560))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file; may be gzipped', required=True)
    parser.add_argument('--output', '-o', help='Output VCF file; may be gzipped', required=True)
    parser.add_argument('--pca', help='Sample and PCA file; may be gzipped', required=True)
    parser.add_argument('--sex', help='Sex (and population file); may be gzipped', required=True)
    parser.add_argument('--consanguineous', help='Consanguinity information; may be gzipped', required=True)
    args = parser.parse_args()
    main(args)