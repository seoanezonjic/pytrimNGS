#! /usr/bin/env python
from pytrimngs.main_modules import *
import argparse

def ls(x): return x.split(',')

def dictio(x): 
    d = {}
    for pair in x.split(';'):
        p = pair.split('=', 1) # split only cuts the first = character
        if len(p) == 2:
            key, val = p
            key = key.strip()
            val = val.strip()
            if ',' in val: val = val.split(',')
        else:
            key = pair[0]
            val = ''            
        d[key] = val
    return d

def pytrimngs(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=" pytrimNGS [OPTIONS]")

    parser.add_argument('-d', '--database',  dest='database', default=None,
        help='Handle databases. "download" for download all databases.')
    parser.add_argument('-t', '--template',  dest='template',
        help='Trimming template to use or path to template')
    parser.add_argument('-Q', '--fastqs',  dest='fastq_files', type=ls,
        help='Input fastq files. Single or paired')
    parser.add_argument('-P', '--parameters', dest='parameters', type=dictio, default={},
        help='Input fastq files. Single or paired')
    parser.add_argument('-w', '--workers',  dest='workers', type=int,
        help='Input fastq files. Single or paired')
    parser.add_argument('-M', '--memory',  dest='memory', 
        help='Input fastq files. Single or paired')
    parser.add_argument('-O', '--output',  dest='output', 
        help='Output folder')
    args = parser.parse_args()
    main_pytrimngs(args)

def pytrimngs_results_parser(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=" pytrimngs_results_parser [OPTIONS]")

    parser.add_argument('-i', '--input_file',  dest='input_file',
        help='pytrimNGS plugin stats directory. This argument must be indicated between quotes')
    args = parser.parse_args()
    main_pytrimngs_results_parser(args)

def get_fastqc_data(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=" get_fastqc_data [OPTIONS]")

    parser.add_argument('-i', '--input_file',  dest='input_file',
        help='fastqc output directory')
    parser.add_argument('-H', '--header',  dest='header', default=False, action='store_true',
        help='Show header')
    parser.add_argument('-m', '--make_mean2count_metrics',  dest='make_mean2count_metrics', default=False, action='store_true',
        help='Make mean across files of the metrics that are absolute counts ("total_sequences')
    parser.add_argument('-T', '--transpose',  dest='transpose', default=False, action='store_true',
        help='Show stat matrix transposed')

    args = parser.parse_args()
    main_get_fastqc_data(args)

def parse_STAR_log(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=" parse_STAR_log [OPTIONS]")

    parser.add_argument('-d', '--data',  dest='data',
        help='Star Log File as input')
    args = parser.parse_args()
    main_parse_STAR_log(args)

def filter_fastq(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description='Filter FASTQ records by sequence length.')
    parser.add_argument('-i', '--input', required=True, dest= 'input',
                        help='Set input file')
    parser.add_argument('-m', '--min_length', type=int, required=True, dest='min_length', 
                        help='Set minimum length size')
    args = parser.parse_args()

    main_filter_fastq(args)

def collapse_bwt(args = None):
    if args == None: args = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Collapse and count identical BWT mappings, updating read names.')
    # Support both -i/--input and a positional argument for compatibility
    parser.add_argument('-i', '--input', help='Input uncollapsed BWT file (tab-delimited)')
    parser.add_argument('input_pos', nargs='?', help=argparse.SUPPRESS)
    args = parser.parse_args()

    main_collapse_bwt(args)

    