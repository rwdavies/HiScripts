#!/usr/bin/env python
# coding:utf-8
import os
import sys
import gzip
import re
import multiprocessing as mp
from calcu_hwe_info_eaf import *


def calculate_info(lines):
    out = []
    for line in lines:
        if re.match('^#', line):
            continue
        tmp = line.rstrip().split('\t')
        gps = array([map(lambda gi:float(gi), g.split(':')[1].split(','))
                     for g in tmp[9:]])

        out.append("\t".join(tmp[0:5]) + "\t" + get_hwe_info_eaf(gps) + "\n")

    return out


def process_wrapper(filename, chunkStart, chunkSize):

    # print 'subprocess pid is %s ' % os.getpid()
    with open(filename) as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
        out = calculate_info(lines)
    return out


def chunkify(filename, size=1024 * 1024):
    fileEnd = os.path.getsize(filename)

    with open(filename, 'r') as f:

        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size, 1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break


def main():

    import argparse
    parser = argparse.ArgumentParser(
        description='Calculate HWE INFO EAF from STITCH output uncompressed vcf with multiprocessing')
    parser.add_argument('input', help='specify the uncompressed vcf')
    parser.add_argument('output', help='output file')
    parser.add_argument('ncores', type=int, help='the number of cores')
    parser.add_argument(
        '-size', type=int, help='the chunk size read from vcf file (default size = 1024 * 1024 * 10 means processing 40 lines per thread for 10k samples)')
    args = parser.parse_args()

    vcf = args.input
    cores = args.ncores
    pool = mp.Pool(cores)
    jobs = []
    size = 1024 * 1024 * 10  # process 40 lines per job for 10k samples
    for chunkStart, chunkSize in chunkify(vcf, size):
        jobs.append(pool.apply_async(
            process_wrapper, (vcf, chunkStart, chunkSize)))

    output = args.output
    o = open(output, 'w')
    for job in jobs:
        o.writelines(job.get())

    pool.close()
    o.close()


if __name__ == '__main__':
    main()
