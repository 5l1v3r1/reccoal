import os
import numpy as np

command = 'make -f Makefile_cluster'
os.system(command)

Llist=[10000] 
Nlist = [3000] #, 10000]
siglist = [0.01, 0.03, 0.1] 
nsamples = 1000
dt=500
script = 'fixed_sigma_burnin.py'
p=0
mem = {300:'1G', 1000:'8G', 3000:'1G', 10000:'8G'}
rt = {300:'0:59:0', 1000:'12:59:0', 3000:'12:59:0', 10000:'12:59:0'}

for L in Llist:
	for N in Nlist:
		for sigma in siglist:
			for rho in np.logspace(np.log10(0.1*sigma),0,5)[:3]:
				arguments = " ".join(map(str,['--pop', N,'--sigma', sigma,
							      '--r', rho, '--L', L]))
				command='qsub -cwd -l h_vmem='+mem[N]+' -l h_rt='+rt[N]+' -p '+str(p)+' submit_script.py '+script+' '+arguments
				os.system(command)
				print command

