import numpy as np
from matplotlib import pyplot as plt
from Bio import Phylo as P
import pickle
import argparse
import glob
from matplotlib import cm

def parse_file_name(fname):
    entries = fname[:-4].split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params

T2={}
T3={}
T4={}
TMRCA={}
diversity={}
SFS={}
file_list = glob.glob('../data/hist_adaptive_*pickle')
for fname in file_list:
	params = parse_file_name(fname)
	N,sigma,rho,L=params['N'], params['Nsigma'],params['Nrho'],params['L']
	with open(fname, 'r') as f:
		data = pickle.load(f)
		label = (N,sigma,rho,L)
	if label in T2:
		T2[label].extend(data['T2'])
		T3[label].extend(data['T3'])
		T4[label].extend(data['T4'])
		TMRCA[label].extend(data['TMRCA'])
                diversity[label].extend(data['selected_diversity'])
                SFS[label][2] += data['SFS'][2]
	else:
		T2[label] = data['T2'].tolist()
		T3[label] = data['T2'].tolist()
		T4[label] = data['T2'].tolist()
		TMRCA[label] = data['TMRCA'].tolist()
                diversity[label] = data['selected_diversity'].tolist()
                SFS[label]=list(data['SFS'])

mT2=np.zeros((len(T2), len(label)+2))
mT3=np.zeros((len(T2), len(label)+2))
mT4=np.zeros((len(T2), len(label)+2))
mTMRCA=np.zeros((len(T2), len(label)+2))
mdiv = np.zeros((len(T2),2))
length = np.zeros(len(T2))
sigma = np.zeros(len(T2))
sigma_dic = {}
for li, label in enumerate(T2.keys()):
	mT2[li,:]=list(label)+[np.mean(T2[label]),np.var(T2[label])]
	mT3[li,:]=list(label)+[np.mean(T3[label]),np.var(T3[label])]
	mT4[li,:]=list(label)+[np.mean(T4[label]),np.var(T4[label])]
	mTMRCA[li,:]=list(label)+[np.mean(TMRCA[label]),np.var(TMRCA[label])]
        mdiv[li,:] = np.mean(diversity[label], axis=0)

	tmproots = np.roots([1, label[2], 0,-label[1]**3])
	sigma[li] = np.min(tmproots[np.where(np.isreal(tmproots)) and np.where(tmproots>0)])
	length[li] = 1.0/(1.0+label[2]/sigma[li])
        sigma_dic[label] = sigma[li]
ii = np.where(mT2[:,1]>0.0)[0]

Nlist = np.unique(mT2[:,0])
siglist = np.unique(mT2[:,1])
rholist = np.unique(mT2[:,2])


######################
# Figure T2
######################
plt.figure()
plt.title('T2: color = LD segment length')
plt.scatter((sigma[ii]/(2*np.log(1+mT2[ii,1]))**0.5), 
            (mT2[ii,-2]/mT2[ii,0]), 
            c=length[ii], s=50)
x=np.logspace(-1.5, 1.5, 20)
plt.plot(x, 1.0/(1+x))
plt.colorbar()
plt.xlabel(r'$N\hat{\sigma}/\sqrt{2\log N\hat{\sigma}}$')
plt.ylabel(r'$T_2/N$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.savefig('../figures/T2_ada.pdf')

######################
# Figure TMRCA
######################
plt.figure()
plt.title(r'$T_{MRCA}$: color = LD segment length')
plt.scatter((sigma[ii]/(2*np.log(1+sigma[ii]))**0.5), 
            (mTMRCA[ii,-2]/mT2[ii,0]), 
            c=(length[ii]), s=50)
x=np.logspace(-1.5, 1.5, 20)
plt.plot(x, 2.0/(1+x))
plt.colorbar()
plt.xlabel(r'$N\hat{\sigma}/\sqrt{2\log N\hat{\sigma}}$')
plt.ylabel(r'$T_{MRCA}/N$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.savefig('../figures/TMRCA_ada.pdf')


###################
# SFS
####################
plt.figure()
plt.title(r'color = $N\sigma$')
for label in SFS:
    plt.plot(np.log(SFS[label][0]/(1-SFS[label][0])), 
             SFS[label][2]/SFS[label][1]/label[0], 
             c=cm.jet((np.log10(sigma_dic[label])-np.log10(sigma.min()))
                      /(np.log10(sigma.max())-np.log10(sigma.min()))))
#             c=cm.jet((np.log10(label[0])-np.log10(Nlist.min()))
#                      /(np.log10(Nlist.max())-np.log10(Nlist.min()))))
#             c=cm.jet((np.log10(label[1])-np.log10(siglist.min()))
#                      /(np.log10(siglist.max())-np.log10(siglist.min()))))

ax=plt.gca()
ax.set_yscale('log')
plt.savefig('../figures/SFS_ada.pdf')         
