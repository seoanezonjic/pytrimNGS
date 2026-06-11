#! /usr/bin/env python
import glob
import sys
import os
import re
import math
import subprocess
from termcolor import colored

############################################################################
## METHODS
############################################################################
def load_template(path):
    parameters = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == '' or re.match(r'^#', line): continue
            variable, value = get_key_val(line)
            parameters[variable] = value
    return parameters

def get_key_val(text):
    pair = text.split('=', 1)
    if len(pair) == 2:
        key, val = pair
        key = key.strip()
        val = val.strip()
        if ',' in val: val = val.split(',')
    else:
        key = pair[0]
        val = ''
    return key, val

def get_db_path(sub_db, db_path):
     return os.path.join(db_path, 'fastas', sub_db, sub_db+'.fasta.gz')

def get_cmd(plugin, out_logs=None):
    cmd = plugin['executor']
    params = plugin['parameters']
    add_opts = plugin.get('additional_parameters')
    if add_opts == None: add_opts = ''

    if params.get('stats'): 
        params['stats'] = os.path.join(out_logs, params['stats'])
    elif params.get('refstats'): 
        params['refstats'] = os.path.join(out_logs, params['refstats'])        
    for param, val in params.items():
        if val !='' and val != None:
            cmd = cmd + f' {param}={val}'
        elif val =='': # To write bool parameters. If False, the pair must be NOT writed
            cmd = cmd + f' {param}'

    cmd = cmd + f" {add_opts}" + f' 2> {os.path.join(out_logs, plugin['output'])}'
    return cmd

def get_full_cmd(plugin_list, plugins, db_path, bb_path_jni, bb_path_current, out_log=None, parms= {}):
    modules = ['PluginReadInputBb']
    modules.extend(plugin_list)
    modules.append('PluginSaveResultsBb')
    cmds = [ ]
    for mod in modules:
        plugin = plugins[mod]
        if mod == 'PluginContaminants':
            for cDB in parms['contaminants_db']:
                if os.path.exists(cDB): #Custom DB
                    dbout = os.path.join(os.path.dirname(cDB), 'index')
                    dbname = os.path.basename(cDB).split('.')[0]
                    index_database(cDB, dbout, bb_path_jni, bb_path_current)
                    plugin['output'] = f"{dbname}_contaminants_filtering_stats_cmd.txt"
                    plugin['parameters']['path'] = dbout
                    plugin['parameters']['refstats'] = f"{dbname}_contaminants_filtering_stats.txt"
                else: # internal DB
                    cDB_path = os.path.join(db_path, 'indices', cDB)
                    plugin['output'] = f"{cDB}_contaminants_filtering_stats_cmd.txt"
                    plugin['parameters']['path'] = cDB_path
                    plugin['parameters']['refstats'] = f"{cDB}_contaminants_filtering_stats.txt"
                add_pluguin_cmd(plugin, mod, cmds, out_log)
        else:
            add_pluguin_cmd(plugin, mod, cmds, out_log)

    pipe_cmd = " | ".join(cmds)
    print(f"\n>> {colored('FULL CMD', 'red')}\n{colored(pipe_cmd, 'blue')}")
    return pipe_cmd

def add_pluguin_cmd(plugin, mod, cmds, out_logs=None):
    cmd = get_cmd(plugin, out_logs=out_logs)
    print(f">> {colored(mod, 'green')}\n{colored(cmd, 'yellow')}")
    cmds.append(cmd)

def execute_cmd(cmd):
    result = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT)
    #result = subprocess.run(cmd.split(' ') capture_output=True, text=True)
    print(f'stdout: {result}')
    #print(f'stdout: {result.stdout.strip()}')
    #print(f'stderr: {result.stderr.strip()}')
    return result

def download_database(urldb, output, current, jni):
    cmd = f"git clone {urldb} {output}"
    print(cmd)
    execute_cmd(cmd)
    index_folder = os.path.join(output, 'indices')
    if not os.path.exists(index_folder): os.mkdir(index_folder)
    for subDBin in glob.glob(os.path.join(output, 'fastas', '*')):
        subDB = os.path.basename(subDBin)
        subDBout = os.path.join(index_folder, subDB)
        index_database(subDBin, subDBout, jni, current)

def index_database(dbin, dbout, jni, current):
    if not os.path.exists(dbout):
        print(f" -- Creating database {dbout} ")
        os.mkdir(dbout)
        cmd = f"java -Djava.library.path={jni} -ea -cp {current} align2.BBSplitter ow=t fastareadlen=500 minhits=1 maxindel=20 qtrim=rl untrim=t trimq=6 t=1 ref={dbin} path={dbout} 2> {os.path.join(dbout, 'log.txt')}" 
        execute_cmd(cmd)
    else:
        print(f" -- The folder database {dbout} exists, skipping indexing. To update, remove the folder before execution")

def get_cpu(plugin_list, contaminants_db, workers):
    cpu_assign = {}
    n_cont_db = 0
    if contaminants_db != None: n_cont_db = len(contaminants_db)
    min_cpu = len(plugin_list) + 2 + n_cont_db -1 # +2 for cpus used in input and output reads. -1 for the cpu already assigned to the PluginContaminant in first iteration
    if workers < min_cpu: print(f'> WARNING! The current configuration needs at least {min_cpu} threads')
    for p in plugin_list: cpu_assign[p] = 1
    free_cpu = workers - min_cpu
    if free_cpu > 0:
        if 'PluginContaminants' in plugin_list:
            cpu4contDB = math.ceil(free_cpu/n_cont_db)
            cpu_assign['PluginContaminants'] += cpu4contDB
        else: #adapters
            adapt_count = 0
            if 'PluginAdapters3' in plugin_list: adapt_count +=1
            if 'PluginAdapters5' in plugin_list: adapt_count +=1
            cpu4adapt = math.ceil(free_cpu/adapt_count)
            cpu_assign['PluginAdapters3'] += cpu4adapt
            cpu_assign['PluginAdapters5'] += cpu4adapt

    print(f"cpu assign: {cpu_assign}") 
    return cpu_assign