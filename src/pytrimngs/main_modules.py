#! /usr/bin/env python
import sys
import os
import re
import time
from importlib.resources import files

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
