import os
import glob
import subprocess as sp

def parse_file_name(fname):
    entries = fname[:-4].split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params


file_list = glob.glob('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_fixed_effect/hist_adaptive_*.pickle')

unique_parameters = set()
for fname in file_list:
	params=parse_file_name(fname)
	unique_parameters.add((params['N'], params['s'], params['r'], params['L']))

for P in unique_parameters:
	sp.call(['python', 'src/analyze_data_fixed_effect.py']+map(str,['--pop', int(P[0]),'--s',P[1],'--r', P[2], '--L', int(P[3])]))


