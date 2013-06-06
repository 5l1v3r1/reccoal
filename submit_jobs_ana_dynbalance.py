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


file_list = glob.glob('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_dynbalance/hist_adaptive_*.pickle')

unique_parameters = set()
for fname in file_list:
	params=parse_file_name(fname)
	unique_parameters.add((params['N'], params['s'], params['U'], params['r'], params['L']))

for P in unique_parameters:
	arguments = " ".join(map(str,['--pop', int(P[0]),'--s',P[1],'--U', P[2],'--r', P[3], '--L', int(P[4])]))
	command='qsub -cwd -l h_rt=0:15:0 submit_script_ana.py '+arguments
#os.system(command)
	sp.call(['python', 'src/analyze_data_dynbalance.py']+map(str,['--pop', int(P[0]),'--s',P[1],'--U', P[2], '--r', P[3], '--L', int(P[4])]))
	print command
	

