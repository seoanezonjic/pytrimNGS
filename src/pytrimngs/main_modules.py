#! /usr/bin/env python
import glob
import sys
import os
import re
import time
from importlib.resources import files

def get_adapters_reads(files_index):
    reads_after_adapters = []
    for file in ["adapters_5_trimming_stats_cmd.txt", "adapters_3_trimming_stats_cmd.txt"]:
        data = files_index.get(file) 
        if data == None:
            reads_after_adapters.append(None)
            continue
        for line in data:
            if line[0] == "Result:": reads_after_adapters.append(int(line[1].split("reads")[0]))

    reads_after_adapters = [x for x in reads_after_adapters if x is not None]
    if len(reads_after_adapters) == 0:
        reads_after_adapters = 'NA'
    else:
        reads_after_adapters = min(reads_after_adapters)
    return reads_after_adapters

def get_ribo_cont_reads(files_index):
    reads_after = 0
    for file in ["ribo_cont_trimming_stats_cmd.txt"]:
        data = files_index.get(file) 
        if data == None:  continue
        for line in data: 
            if line[0] == "Result:": reads_after += int(line[1].split("reads")[0])
    return reads_after

def get_contaminant_reads(files_index):
    all_contaminants = 0
    for file in ["contaminants_contaminants_filtering_stats.txt"]:
        data = files_index.get(file)
        if data != None:
            for contaminant_data in data:
                if contaminant_data[0] == "name": continue
                all_contaminants += int(contaminant_data[7])
    return all_contaminants


def main_pytrimngs_results_parser(opts):
    args = vars(opts)
    files_index = {}
    for file_path in glob.glob(os.path.join(args['input_file'], '*.txt')):
        filename = os.path.basename(file_path)
        file_data = []
        with open(file_path) as f:
            for line in f: 
                line = line.rstrip()
                line = re.sub("#","", line)
                line = re.sub(r' +', "", line)
                file_data.append(line.split("\t"))
            files_index[filename] = file_data
    
    stbb_metrics = {}
    stbb_metrics["adapter_filter_passed"] = get_adapters_reads(files_index)
    stbb_metrics["contaminants_reads"] = get_contaminant_reads(files_index)
    stbb_metrics["ribo_contaminants_reads"] = get_ribo_cont_reads(files_index)
    print_metrics(stbb_metrics)

def print_metrics(stbb_metrics):
    for metric_name, metric in stbb_metrics.items():
        print(f"{metric_name}\t{metric}")

def main_get_fastqc_data(args):
    import zipfile
    import io
    from pytrimngs.fastqc_parser import FastQC_Parser
    options = vars(args)
    all_stats = []
    header = ['total_sequences', 'read_max_length', 'read_min_length', '%gc', 'mean_qual_per_base', 
              'min_qual_per_base_in_lower_quartile', 'min_qual_per_base_in_10th_decile', 
              'weigthed_qual_per_sequence', 'mean_indeterminations_per_base', 
              'weigthed_read_length', 'sequence_length_distribution']
    for file in glob.glob(options['input_file']):
        with zipfile.ZipFile(file) as zf:
            zip_base_folder = zf.namelist()[0] 
            text_file_handler = zf.open(name = os.path.join(zip_base_folder, "fastqc_data.txt"), mode = 'r')
            text_file = io.TextIOWrapper(text_file_handler, encoding='utf-8', newline='').read()
            modules = FastQC_Parser.parse_fastqc_data(text_file)
            stats = [
                modules['general_stats']['Total Sequences'],
                modules['general_stats']['Read_max_length'],
                modules['general_stats']['Read_min_length'],
                modules['general_stats']['%GC'],
                FastQC_Parser.get_mean(modules['quality_per_base'], 2),
                FastQC_Parser.get_min(modules['quality_per_base'], 3),
                FastQC_Parser.get_min(modules['quality_per_base'], 5),
                FastQC_Parser.get_weighted_mean(modules['quality_per_sequence'], 0,1),
                FastQC_Parser.get_mean(modules['indeterminations_per_base'], 1),
                FastQC_Parser.get_weighted_mean_with_intervals(modules['sequence_length_distribution'], 0,1),
                FastQC_Parser.parse_distributions(modules['sequence_length_distribution'])
            ]
            all_stats.append(stats)

    n_samples = len(all_stats)
    n_parameters = len(header)
    means = []
    for parameter_index in range(0, n_parameters):
        if header[parameter_index] == 'sequence_length_distribution':
            all_distributions = [ all_stats[s][parameter_index] for s in range(0, n_samples) ]
            means.append(";".join(all_distributions))
        else:
            summ = 0
            for sample_index in range(0, n_samples): summ += all_stats[sample_index][parameter_index]
            if not options['make_mean2count_metrics'] and header[parameter_index] == 'total_sequences':
                means.append(summ)
            else:
                means.append(summ/n_samples)
    
    means_string = [ str(m) for m in means ]
    if not options['transpose']:
        if options['header']: print("\t".join(header))
        print("\t".join(means_string))
    else:
        for i, mean in enumerate(means_string):
            record = []
            if options['header']: record.append(header[i])
            record.append(mean)
            print("\t".join(record))

def main_parse_STAR_log(args):
    opts = vars(args)
    with open(opts['data']) as f:
        for line in f:
            if '|' in line:
                metric, value = line.strip().split("\t")
                metric = re.sub('\\|', '', metric)
                metric = metric.strip()
                metric = re.sub(' ', '_', metric)
                print(f"{metric}\t{value}")

def main_filter_fastq(args):
    opts = vars(args)

    entry = []
    attributes_parsed = 0
    with open(opts['input'], 'r') as fh:
        for line in fh:
            line = line.rstrip()
            attributes_parsed += 1
            entry.append(line)
            if attributes_parsed == 4:
                if len(entry[1]) >= opts['min_length']: print("\n".join(entry))
                entry = []
                attributes_parsed = 0

READ_NAME = 0
STRAND = 1
CHR = 2
COORD = 3
SEQ = 4

def main_collapse_bwt(opts):
    args = vars(opts)

    input_bwt = args['input'] if args['input'] is not None else args['input_pos']
    if not input_bwt or not os.path.isfile(input_bwt) or os.stat(input_bwt).st_size == 0:
        sys.stderr.write("No input file set\n\n\tUSAGE:\tcollapse_bwt uncollapsed.bwt > collapsed.bwt\n")
        sys.exit(1)

    # Use streaming indexing approach
    indexed_bwt = load_and_index_bwt(input_bwt)
    collapsed_bwt = build_collapsed_bwt(indexed_bwt)
    print_output(collapsed_bwt)

def load_and_index_bwt(input_file):
    """
    Read file line-by-line and create an ordered index:
    key -> [first_mapping, count]
    The insertion order is preserved .
    """
    indexed_bwt = {}
    with open(input_file, 'r') as fh:
        for line in fh:
            mapping = line.rstrip('\n').split('\t')
            key = f"{mapping[CHR]}:{mapping[STRAND]}:{mapping[COORD]}:{mapping[SEQ]}"
            if key not in indexed_bwt:
                indexed_bwt[key] = [mapping, 1]
            else:
                indexed_bwt[key][1] += 1
    return indexed_bwt

def build_collapsed_bwt(indexed_bwt):
    """
    Build collapsed entries. For each indexed key (in insertion order),
    update stacked_mappings and change the read name to:
      seq_{stacked_mappings}_x{collapsed_mappings}
    """
    collapsed_bwt = []
    stacked_mappings = 0
    for _id, attributes in indexed_bwt.items():
        bwt_entry, collapsed_mappings = attributes
        stacked_mappings += collapsed_mappings
        collapsed_bwt.append(change_read_name(stacked_mappings, collapsed_mappings, bwt_entry))
    return collapsed_bwt

def change_read_name(stacked_mappings, collapsed_mappings, bwt_entry):
    """Replace READ_NAME field with the collapsed name and return the entry."""
    bwt_entry[READ_NAME] = f"seq_{stacked_mappings}_x{collapsed_mappings}"
    return bwt_entry

def print_output(bwt):
    """Print each bwt entry as a tab-separated line."""
    for bwt_entry in bwt: print("\t".join(bwt_entry))    

def mask_fasta(mask_regs, path):
    with open(path) as f:
        s_id = None
        seq = ''
        for line in f:
            line = line.rstrip()
            if re.search('^>', line):
                line = line.replace('>', '')
                regs = mask_regs.get(s_id)
                if regs != None: seq = mask_seq(seq, regs)
                print(f">{s_id}\n{seq}")
                s_id = line
                seq = ''
            else:
                seq = seq + line
    regs = mask_regs.get(s_id)
    if regs != None: seq = mask_seq(seq, regs)
    print(f">{s_id}\n{seq}")

#bed file formated taking into account info in https://www.ensembl.org/info/genome/genebuild/human_PARS.html
def load_bed(path):
    regs = {}
    with open(path) as f: 
        for line in f: 
            chrm_id, start, stop = line.rstrip().split("\t")
            chrm = regs.get(chrm_id)
            if chrm == None:
                chrm = []
                regs[chrm_id] = chrm
            chrm.append([int(start) -1, int(stop) -1])
    return regs

def mask_seq(seq, regs):
    new_seq = ''
    current_coord = 0
    for start, stop in regs:#slicing gets +1 position
        new_seq = new_seq + seq[current_coord:start] + 'N' * (stop +1 - start)
        current_coord = stop + 1
    new_seq = new_seq + seq[current_coord:]
    return new_seq

def main_maskFasta(opts):
    args = vars(opts)   
    seq_ids = load_bed(args['bed_file'])
    mask_fasta(seq_ids, args['regions'])

def filter_fasta(ids, path):
    with open(path) as f:
        s_id = None
        seq = ''
        for line in f:
            line = line.rstrip()
            if re.search('^>', line):
                line = line.replace('>', '')
                if s_id in ids: print(f">{s_id}\n{seq}")
                s_id = line
                seq = ''
            else:
                seq = seq + line
    if s_id in ids: print(f">{s_id}\n{seq}")

def load_ids(path):
    ids = []
    with open(path) as f:
        for line in f: ids.append(line.rstrip())
    return ids

def main_lista2fasta(opts):
    args = vars(opts)   
    seq_ids = load_ids(args['seq_ids'])
    filter_fasta(seq_ids, args['fasta'])

def main_get_too_short(opts):
    import pysam
    args = vars(opts) 
    bamfile = pysam.AlignmentFile(args['input'], "rb")
    outfile = pysam.AlignmentFile(args['output'] + ".bam", "w", template=bamfile)
    outFasta = open(args['output'] + ".fasta",'w')
    unaligned_reads = 0
    for read in bamfile:#.fetch('chr1'):
        if read.is_secondary or read.is_supplementary: continue
        read_len = read.infer_read_length()
        cigar_nucleotides, cigar_blocks = read.get_cigar_stats()
        #print(f"{cigar_nucleotides[0]} {cigar_nucleotides[4]} {read_len}")
        soft_clipping_rate = cigar_nucleotides[4]/read_len
        if soft_clipping_rate >= args['min_soft_clip_rate']:
            outfile.write(read)
            pair_read = 1
            if read.is_read2: pair_read = 2
            outFasta.write(f">{read.query_name}_{pair_read}\n{read.query_sequence}\n")
            unaligned_reads += 1

    print(f"unaligned_reads: {unaligned_reads}")
    outFasta.close()    

def main_pytrimngs(opts):
    import pytrimngs # For external_data
    from pytrimngs.pytrimngs import load_template, get_key_val, get_db_path, get_cmd, get_full_cmd, execute_cmd, download_database, index_database, get_cpu 

    args = vars(opts)

    #If there is a BBDB environment var (databases location), then use it
    if os.environ.get('BBDB') != None:
	    db_path = os.environ['BBDB']
    else: # otherwise use SEQTRIM_PATH + DB
	    db_path = str(files('pytrimngs').joinpath('DB'))

    #First set a BBtools path, then checks if BBtools is properly installed. If there is a BBTOOLS_PATH environment var, then use it
    if os.environ.get('BBTOOLS_PATH') != None:
	    bb_path = os.environ['BBTOOLS_PATH']
    else: # otherwise use the result of which bbmap.sh
        bb_path = os.path.dirname(execute_cmd("which bbmap.sh"))
    bb_path_current = os.path.join(bb_path, 'current')
    bb_path_jni = os.path.join(bb_path, 'jni')

    if args['database'] == 'download':
        if not os.path.exists(db_path): os.mkdir(db_path)
        download_database('https://github.com/seoanezonjic/pytrimngs-databases.git', db_path, bb_path_current, bb_path_jni)
        quit()

    #Load template
    if os.path.exists(args['template']):
        template_path = args['template']
    else:
        template_path = str(files('pytrimngs.templates').joinpath(args['template']))
    print(template_path)

    # input/output
    tmp_output = os.path.join(args['output'], 'output_files_tmp')
    out_logs = os.path.join(tmp_output, 'plugins_logs')
    if not os.path.exists(tmp_output):
        os.mkdir(tmp_output)
        os.mkdir(out_logs)
    else:
        sys.exit("output_files_tmp folder already exists so there is a broken execution") 

    final_output = os.path.join(args['output'], 'output_files')
    infastq1 = args['fastq_files'][0]
    if len(args['fastq_files']) == 2:
        infastq2 = args['fastq_files'][1]
        outfastq1 = os.path.join(tmp_output, 'paired_1.fastq.gz') 
        outfastq2 = os.path.join(tmp_output, 'paired_2.fastq.gz')
        interleaved = 't' 
    else:
        infastq2 = False
        outfastq1 = os.path.join(tmp_output, 'single_.fastq.gz') 
        outfastq2 = False 
        interleaved = 'f' 


    parms = {
        'adapters_3_max_mismatches': 1, # PluginAdapters3
        'adapters_5_max_mismatches': 1,  # PluginAdapters5
    }
    parms.update(load_template(template_path))
    parms.update(args['parameters'])
    if parms.get('plugin_list') != None and not isinstance(parms['plugin_list'], list): parms['plugin_list'] = [parms['plugin_list']]
    if parms.get('contaminants_db') != None and not isinstance(parms['contaminants_db'], list): parms['contaminants_db'] = [parms['contaminants_db']]
    cpu_assignation = get_cpu(parms['plugin_list'], parms.get('contaminants_db'), args['workers'])
    # Check adapters DB
    adapters_db = parms.get('adapters_db')
    if adapters_db == None:
        adapters_db = get_db_path('adapters', db_path)
    else:
        dbout = os.path.join(os.path.dirname(adapters_db), 'index')
        index_database(adapters_db, dbout, bb_path_jni, bb_path_current)
    cont_rib_base = os.path.join(db_path, 'fastas', 'cont_ribosome')
    cont_rib_db = f"{os.path.join(cont_rib_base, 'rrna_lsu90.fasta.gz')},{os.path.join(cont_rib_base, 'rrna_ssu90.fasta.gz')}"
    plugins = {
        'PluginReadInputBb': {
            'executor': f'java -ea -cp {bb_path_current} jgi.ReformatReads t=1 -Xmx{args['memory']}m',
            'parameters': {'in': infastq1, 'int': 'f', 'out': 'stdout.fastq'}, # 'in2': infastq2,
            'output': 'input_stats.txt'
            },
        'PluginAdapters3': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginAdapters3')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                    'stats': 'adapters_3_trimming_stats.txt', 'k': 15, 'mink':8, 'hdist': parms['adapters_3_max_mismatches'], 'ktrim': 'r', 'ref': adapters_db,  'tbo': '', 'tpe': ''},
            'output': 'adapters_3_trimming_stats_cmd.txt'
        },
        'PluginAdapters5': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginAdapters5')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                    'stats': 'adapters_5_trimming_stats.txt', 'k': 21, 'mink':11, 'hdist': parms['adapters_5_max_mismatches'], 'ktrim': 'l', 'ref': adapters_db,  'tbo': '', 'tpe': ''},
            'output': 'adapters_5_trimming_stats_cmd.txt'
        },
        'PluginPolyAt': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginPolyAt')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                       'trimpolya': 9 },
            'output': 'polyat_trimming_stats.txt'
        },
        'PluginContaminants': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} align2.BBSplitter t={cpu_assignation.get('PluginContaminants')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved, 
                       'ow': 't', 'fastareadlen': 500, 'minhits': 1, 'maxindel': 20, 'qtrim': 'rl',
                        'untrim': 't', 'trimq': 6, 'minratio': 0.56}, # path=/indices/contaminants refstats=/output_files_tmp
            'additional_parameters': parms.get('contaminants_aditional_params'), 
            'output': None
        },
        'PluginRiboContaminants': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginRiboContaminants')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved, 
                       'stats': 'ribo_cont_trimming_stats.txt', 'ref': cont_rib_db, 'k':31 },
            'output': 'ribo_cont_trimming_stats_cmd.txt'
        },        
        'PluginQuality':{
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginQuality')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                       'trimq': parms['quality_threshold'], 'qtrim': 'rl'},
            'additional_parameters': parms.get('quality_aditional_params'), 
            'output': 'quality_trimming_stats.txt'
        },
        'PluginLowComplexity':{
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={cpu_assignation.get('PluginLowComplexity')} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                       'entropy': 0.01, 'entropywindow': 50, 'minlength': 50},
            'output': 'low_complexity_stats.txt'
        },
        'PluginSaveResultsBb': {
            'executor': f'java -ea -cp {bb_path_current} jgi.ReformatReads t=2 -Xmx{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': interleaved,
                           'out': outfastq1, 'minlength': parms['minlength']}, #'out2': outfastq2,
            'output': 'output_stats.txt'
        }
    }
    # USER MODIFICATIONS
    # In single end, disable IN/OUT of second pair
    if infastq2 != False: plugins['PluginReadInputBb']['parameters']['in2'] = infastq2
    if outfastq2 != False: plugins['PluginSaveResultsBb']['parameters']['out2'] = outfastq2

    pipe_cmd = get_full_cmd(parms['plugin_list'], plugins, db_path, bb_path_jni, bb_path_current, out_log=out_logs, parms=parms)
    execute_cmd(pipe_cmd)
    os.rename(tmp_output, final_output)
