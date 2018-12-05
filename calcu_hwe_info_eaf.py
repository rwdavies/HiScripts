#!/usr/bin/env python
import gzip
import sys
import re
import math
from numpy import *


def get_max_gen_rapid(gps):

    (nrow, ncol) = gps.shape
    d0 = full((nrow, 2), False)
    d1 = full((nrow, 1), True)
    d = hstack((d1, d0))
    x = gps.copy()
    y = x[:, 0]
    for i in [1, 2]:
        w = x[:, i] > y
        d[w, i] = True
        d[w, i - 1] = False
        d[w, 0] = False
        y[w] = x[w, i]

    return(d)


def calculate_hwe_p(hweCount):
    # x = array([100,50,10])
    x = sum(hweCount, axis=0)
    if x[2] > x[0]:
        x = x[[2, 1, 0]]
    nAA = x[0]
    nAB = x[1]
    nBB = x[2]
    nA = nAA * 2 + nAB
    nB = nBB * 2 + nAB
    n = nAA + nAB + nBB  # number of (diploid) individuals

    min_het = 0
    max_het = int(nAB + 2 * min(nAA, nBB))

    # maximum-ish value
    mid = int(floor((nA * nB) / (nA + nB)))
    if math.isnan(mid):
        mid = 0
    # make odd if it needs to be
    if nA % 2 != mid % 2:
        mid = mid + 1

    # determine a mid point
    probs = zeros(shape=(max_het + 1, 1), dtype=float)
    probs[mid] = 1
    # use recurrence relation - going down
    n_het = mid
    n_hom_alt = (nBB * 2 + nAB - n_het) / 2
    n_hom_ref = n - n_het - n_hom_alt
    if mid - 2 >= min_het:
        for het in range(mid - 2, min_het - 2, -2):
            probs[het] = probs[het + 2] * n_het * \
                (n_het - 1) / (4 * (n_hom_ref + 1) * (n_hom_alt + 1))
            n_het = n_het - 2
            n_hom_ref = n_hom_ref + 1
            n_hom_alt = n_hom_alt + 1

    # use recurrence ralationship - going up
    n_het = mid
    n_hom_alt = (nBB * 2 + nAB - n_het) / 2
    n_hom_ref = n - n_het - n_hom_alt
    if mid + 2 <= max_het:
        for het in range(mid + 2, max_het + 2, 2):
            probs[het] = probs[het - 2] * \
                (4 * n_hom_ref * n_hom_alt) / ((n_het + 2) * (n_het + 1))
            n_het = n_het + 2
            n_hom_ref = n_hom_ref - 1
            n_hom_alt = n_hom_alt - 1

    all_probs = probs / sum(probs, axis=0)
    p2 = sum(all_probs[all_probs <= all_probs[int(nAB)]])
    p = min(1, p2)

    return(p)


def get_hwe_info_eaf(gps):

    seterr(divide='ignore', invalid='ignore', over='ignore')
    # info count
    N = gps.shape[0]
    eij = gps[:, 1] + gps[:, 2] * 2
    fij = gps[:, 1] + gps[:, 2] * 4
    infoCount = zeros(shape=(N, 2))
    infoCount[:, 0] = eij + infoCount[:, 0]
    infoCount[:, 1] = fij - eij * eij + infoCount[:, 1]
    allinfoCount = sum(infoCount, axis=0)
    thetaHat = allinfoCount[0] / 2 / N
    denom = 2 * N * thetaHat * (1 - thetaHat)
    info = around(1 - allinfoCount[1] / denom, 5)
    if around(thetaHat, 2) == 0 or around(thetaHat, 2) == 1:
        info = 1
    elif info < 0:
        info = 0
    else:
        info = info

    # af count
    afCount = zeros(shape=(N, 1))
    afCount[:, 0] = afCount[:, 0] + eij / 2
    allafCount = sum(afCount, axis=0)
    eaf = around(allafCount[0] / N, 5)

    # hwe count
    hweCount = zeros(shape=(N, 3), dtype=float)
    w = get_max_gen_rapid(gps)
    hweCount[w] = hweCount[w] + 1
    hwe_p = calculate_hwe_p(hweCount)

    # output
    out = "EAF=" + str(eaf) + ";INFO_SCORE=" + str(info) + \
        ";HWE=" + str(hwe_p)
    return(out)


def main():
    vcf = sys.argv[1]
    output = sys.argv[2]
    o = open(output, 'w')
    with gzip.open(vcf, 'rb') as f:
        for line in f:
            if re.match('^##', line):
                continue
            tmp = line.rstrip().split('\t')

            if re.match('^#CHROM', line):
                out = "\t".join(['CHROM', 'POS', 'REF', 'ALT'] + tmp[9:])
            else:
                gps = array(
                    [map(lambda gi:float(gi), g.split(':')[1].split(',')) for g in tmp[9:]])
                out = "\t".join(tmp[0:2] + tmp[3:5]) + \
                    "\t" + get_hwe_info_eaf(gps)

            o.write(out + "\n")
    o.close()


if __name__ == '__main__':
    main()
