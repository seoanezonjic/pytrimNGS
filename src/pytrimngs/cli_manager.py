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

