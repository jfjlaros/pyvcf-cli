#!/usr/bin/env python

from __future__ import division

import argparse

import vcf
import vcf.utils as vcfutils

def diff(input_handles, output_handle, precision=10):
    """
    """
    first_vcf = vcf.Reader(input_handles[0])
    second_vcf = vcf.Reader(input_handles[1])
    symmetric_difference = 0
    total = 0

    walker = vcfutils.walk_together(first_vcf, second_vcf)
    for first_record, second_record in walker:
        if first_record and second_record and not (first_record.is_indel or
                second_record.is_indel):
            if (first_record.alleles[1].sequence !=
                    second_record.alleles[1].sequence):
                symmetric_difference += 1
            total +=1
        #if
    #for
    output_handle.write('{value:.{precision}f}\n'.format(
        value=symmetric_difference / total, precision=precision))
#diff

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_handles', metavar='INPUT',
        type=argparse.FileType('r'), nargs=2, help='VCF file')
    parser.add_argument('-o', dest='output_handle',
        type=argparse.FileType('w'), default='-', help='output file')
    parser.add_argument('-p', dest='precision', type=int, default=10,
        help='precision (%(type)s default=%(default)s)')
    parser.set_defaults(func=diff)

    try:
        arguments = parser.parse_args()
    except IOError, error:
        parser.error(error)

    try:
        arguments.func(**{k: v for k, v in vars(arguments).items()
            if k not in ('func', 'subcommand')})
    except ValueError, error:
        parser.error(error)
#main

if __name__ == '__main__':
    main()
