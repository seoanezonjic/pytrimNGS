#! /usr/bin/env python
from pytrimngs.main_modules import *
import argparse

def ls(x): return x.split(',')

def dictio(x): 
    d = {}
    for pair in x.split(';'):
        k, v = pair.split('=')
        d[k] = v
    return d

def pytrimngs(args = None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=" pytrimNGS [OPTIONS]")

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


