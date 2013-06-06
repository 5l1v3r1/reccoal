import os
import glob
import subprocess as sp

def parse_file_name(fname):
    entries = fname.split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params

#data_dir = 'data_dynbalance'
data_dir = 'data_fixed_sigma'
#data_dir = 'data_fixed_effect'
analysis_script = 'src/consolidate_data.py'

file_list = glob.glob('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/'+data_dir+'/hist_adaptive_*.pickle')

unique_parameters = set()
for fname in file_list:
    params=parse_file_name(fname[:-7])
    arg_list = ['--pop', int(params['N']),'--r', params['r'], '--L', int(params['L'])]
    if 'sigma' in params: arg_list.extend(['--sigma', params['sigma']])
    elif 's' in params: arg_list.extend(['--s', params['s']])
    if 'U' in params: arg_list.extend(['--U', params['U']])
    unique_parameters.add(tuple(arg_list))

for arg_list in unique_parameters:
    arguments = " ".join(map(str,arg_list))
    print arguments
    sp.call(['python', analysis_script]+map(str,arg_list))
	

