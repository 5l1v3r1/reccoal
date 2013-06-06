import os
import glob
import numpy as np

def parse_file_name(fname):
    entries = fname[:-4].split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params


file_list = glob.glob('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_new/hist_adaptive_*.pickle')
target_file_size = 9e5
nsamples = 1000
dt=500
file_sizes = {}
mem = {300:'1G', 1000:'8G', 3000:'1G', 10000:'8G'}
rt = {300:'0:59:0', 1000:'0:59:0', 3000:'0:59:0', 10000:'3:59:0'}

for fname in file_list:
	params=parse_file_name(fname)
	if 'Nr' in params: Nr=params['Nr']
	else: Nr=params['Nrho']
	label = (params['N'], params['Nsigma'], Nr, params['L'])
	st = os.stat(fname)
	if (label in file_sizes):
		file_sizes[label].append(st.st_size)
	else:
		file_sizes[label]=[st.st_size]


for label in file_sizes:
	missing_size = target_file_size - np.sum(file_sizes[label])
	run_size = np.mean(file_sizes[label])
	nruns = int(missing_size/run_size)
	N=int(label[0])
	sigma = label[1]/N
	r=label[2]/N
	L=int(label[3])
	for ri in range(300,300+nruns):
		arguments = " ".join(map(str,['--pop', N,'--sigma', sigma,
								  '--r', r, '--L', L, 
								  '--dt', dt, '--n', 100, 
								  '--id', ri]))
		command='qsub -cwd -l h_vmem='+mem[N]+' -l h_rt='+rt[N]+' submit_script_adaptive.py '+arguments
		os.system(command)
		print command
	
