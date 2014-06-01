#!/usr/bin/env python

"""
VCF manipulation toolkit.


Copyright (c) 2014 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2014 Jeroen F.J. Laros <J.F.J.Laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""

from __future__ import division

import argparse
import os

import wiggelen

import vcf
import vcf.utils as vcfutils

__version_info__ = ('0', '0', '1')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://github.com/jfjlaros/PyVCF'

usage = __doc__.split("\n\n\n")

def doc_split(func):
    return func.__doc__.split("\n\n")[0]

def version(name):
    return "%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s" % (name,
        __version__, __author__, __contact__, __homepage__)

def diff(input_handles, output_handle, precision=10):
    """
    Calculate the Jaccard distance between two VCF files.


    :arg input_handles: List of two open readable handles to VCF files.
    :type input_handles: list(stream)
    :arg output_handle: An open writable handle.
    :type output_handle: stream
    :arg precision: Number of decimals in the output.
    :type precision: int
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

def vcf2wig(input_handle, output_handle, field="AF", prefix="",
        function="min"):
    """
    Convert a VCF file to a wiggle track.

    :arg input_handle: Open readable handle to a VCF file.
    :type input_handle: stream
    :arg output_handle: Open writable handle to a wiggle track.
    :type output_handle: stream
    """
    def unpack(list_or_value):
        if type(list_or_value) == list:
            return func(list_or_value)
        return list_or_value

    reader = vcf.Reader(input_handle)
    func = eval("lambda l: {}(l)".format(function))

    wiggelen.write(map(lambda x: ("{}{}".format(prefix, x.CHROM), x.POS,
        x.INFO.has_key(field) and unpack(x.INFO[field]) or 0.0), reader),
        track=output_handle,
        name=os.path.splitext(os.path.basename(output_handle.name))[0])
#vcf2wig

def main():
    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="input file")

    pair_in_parser = argparse.ArgumentParser(add_help=False)
    pair_in_parser.add_argument("input_handles", metavar="INPUT", nargs=2,
        type=argparse.FileType('r'), help="pair of input files")

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("output_handle", metavar="OUTPUT",
        type=argparse.FileType('w'), default='-', help="output file")

    prefix_parser = argparse.ArgumentParser(add_help=False)
    prefix_parser.add_argument("-x", dest="prefix", type=str, default="",
        help='prefix for chromosome names (%(type)s default="%(default)s")')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers()

    parser_diff = subparsers.add_parser("diff", parents=[pair_in_parser,
        output_parser], description=doc_split(diff))
    parser_diff.add_argument('-p', dest='precision', type=int, default=10,
        help='precision (%(type)s default=%(default)s)')
    parser_diff.set_defaults(func=diff)

    parser_vcf2wig = subparsers.add_parser("vcf2wig", parents=[input_parser,
        output_parser, prefix_parser], description=doc_split(vcf2wig))
    parser_vcf2wig.add_argument('-f', dest='field', type=str, default="AF",
        help='INFO field to convert (%(type)s default=%(default)s)')
    parser_vcf2wig.set_defaults(func=vcf2wig)

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
