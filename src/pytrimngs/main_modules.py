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
              'weigthed_read_length sequence_length_distribution']
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
    for parameter_index in range(0, n_parameters - 1):
        if header[parameter_index] == 'sequence_length_distribution':
            all_distributions = [ all_stats[s][parameter_index] for s in range(0, n_samples - 1) ]
            means.append(all_distributions.join(";"))
        else:
            summ = 0
            for sample_index in range(0, n_samples -1): summ += all_stats[sample_index][parameter_index]
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

def main_pytrimngs(opts):
    import pytrimngs # For external_data
    from pytrimngs.pytrimngs import load_template, get_key_val, get_db_path, get_cmd, get_full_cmd, execute_cmd 

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
    else:
        infastq2 = False
        outfastq1 = os.path.join(tmp_output, 'single_.fastq.gz') 
        outfastq2 = False 


    parms = load_template(template_path)
    parms.update(args['parameters'])
    plugins = {
        'PluginReadInputBb': {
            'executor': f'java -ea -cp {bb_path_current} jgi.ReformatReads t={args['workers']} -Xmx{args['memory']}m',
            'parameters': {'in': infastq1, 'in2': infastq2, 'int': 'f', 'out': 'stdout.fastq'},
            'output': 'input_stats.txt'
            },
        'PluginAdapters3': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={args['workers']} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                    'stats': 'adapters_3_trimming_stats.txt', 'k': 15, 'mink':8, 'hdist':1, 'ktrim': 'r', 'ref': get_db_path('adapters', db_path),  'tbo': '', 'tpe': ''},
            'output': 'adapters_3_trimming_stats_cmd.txt'
        },
        'PluginAdapters5': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={args['workers']} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                    'stats': 'adapters_5_trimming_stats.txt', 'k': 21, 'mink':11, 'hdist':1, 'ktrim': 'l', 'ref': get_db_path('adapters', db_path),  'tbo': '', 'tpe': ''},
            'output': 'adapters_5_trimming_stats_cmd.txt'
        },
        'PluginPolyAt': {
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={args['workers']} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                       'trimpolya': 9 },
            'output': 'polyat_trimming_stats.txt'
        },
        'PluginQuality':{
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={args['workers']} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                       'trimq': parms['quality_threshold'], 'qtrim': 'rl'}, 
            'output': 'quality_trimming_stats.txt'
        },
        'PluginLowComplexity':{
            'executor': f'java -Djava.library.path={bb_path_jni} -ea -cp {bb_path_current} jgi.BBDuk t={args['workers']} -Xmx{args['memory']}m -Xms{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                       'entropy': 0.01, 'entropywindow': 50, 'minlength': 50},
            'output': 'low_complexity_stats.txt'
        },
        'PluginSaveResultsBb': {
            'executor': f'java -ea -cp {bb_path_current} jgi.ReformatReads t={args['workers']} -Xmx{args['memory']}m',
            'parameters': {'in': 'stdin.fastq', 'out': 'stdout.fastq', 'int': 't',
                           'out': outfastq1, 'out2': outfastq2, 'minlength': parms['minlength']},
            'output': 'output_stats.txt'
        }
    }

    pipe_cmd = get_full_cmd(parms['plugin_list'], plugins, out_log=out_logs)
    execute_cmd(pipe_cmd)
    os.rename(tmp_output, final_output)
