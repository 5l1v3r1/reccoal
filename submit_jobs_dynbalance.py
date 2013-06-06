import os
import numpy as np

command = 'make -f Makefile_cluster'
os.system(command)

Llist=[3000, 10000] 
Nlist = [1000,3000, 10000]
slist = [0.001, 0.003, 0.01] 
Ulist = [1.0,3.0, 10.0, 30.0]
nsamples = 1000
dt=500
script = 'dynbalance_burnin.py'
p=0
mem = {300:'1G', 1000:'8G', 3000:'1G', 10000:'8G'}
rt = {300:'0:59:0', 1000:'12:59:0', 3000:'12:59:0', 10000:'12:59:0'}

for L in Llist:
	for N in Nlist:
		for s in slist:
			for U in Ulist:
				for rho in np.logspace(np.log10(s),0,5):
					arguments = " ".join(map(str,['--pop', N,'--s', s, '--U', U, '--r', rho, '--L', L]))
					command='qsub -cwd -l h_vmem='+mem[N]+' -l h_rt='+rt[N]+' -p '+str(p)+' submit_script.py '+script+' '+arguments
					os.system(command)
					print command

