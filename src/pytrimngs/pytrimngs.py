#! /usr/bin/env python
import glob
import sys
import os
import re
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
    pair = text.split('=')
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
    for param, val in params.items():
        if val !='' and val != None:
            cmd = cmd + f' {param}={val}'
        elif val =='': # To write bool parameters. If False, the pair must be NOT writed
            cmd = cmd + f' {param}'

    cmd = cmd + f" {add_opts}" + f' 2> {os.path.join(out_logs, plugin['output'])}'
    return cmd

def get_full_cmd(plugin_list, plugins, out_log=None):
    modules = ['PluginReadInputBb']
    modules.extend(plugin_list)
    modules.append('PluginSaveResultsBb')

    cmds = [ ]
    for mod in modules:
        cmd = get_cmd(plugins[mod], out_logs=out_log)
        print(f">> {colored(mod, 'green')}\n{colored(cmd, 'yellow')}")
        cmds.append(cmd)

    pipe_cmd = " | ".join(cmds)
    print(f"\n>> {colored('FULL CMD', 'red')}\n{colored(pipe_cmd, 'blue')}")
    return pipe_cmd

def execute_cmd(cmd):
    result = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT)
    #result = subprocess.run(cmd.split(' ') capture_output=True, text=True)
    print(f'stdout: {result}')
    #print(f'stdout: {result.stdout.strip()}')
    #print(f'stderr: {result.stderr.strip()}')
    return result

def download_database(urldb, output, current, jni):
    cmd = f"git clone {urldb} {output}"
    index_folder = os.path.join(output, 'indices')
    if not os.path.exists(index_folder): os.mkdir(index_folder)
    execute_cmd(cmd)
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