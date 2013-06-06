import numpy as np
from Bio import Phylo as P
import pickle
import argparse
import glob

def rsq(nu12, nu1, nu2):
    return (nu12-nu1*nu2)**2/(nu1*(1-nu1)*nu2*(1-nu2))


#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a purifying selection with recombination")
parser.add_argument('--pop', default=10000, type=int, help='Population size (N)')
parser.add_argument('--L', default=10000, type=int, help='number of loci (L)')
parser.add_argument('--r', default=0.0, type=float, help='crossover rate per genome')
parser.add_argument('--sigma', default=-1, type=float, help='fitness std dev sigma')
parser.add_argument('--s', default=-1, type=float, help='fitness std dev sigma')
parser.add_argument('--U', default=-1, type=float, help='mutation rate per genome')
params=parser.parse_args()

def parse_file_name(fname):
    entries = fname.split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params

LD_func=rsq
N=params.pop
L=params.L
sigma=params.sigma
r=params.r
U=params.U

label = (N,Nsigma,Nr,L)

T2=[]
T3=[]
T4=[]
TMRCA=[]
diversity=[]
pair_freq=[]
fmask = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_new/hist_adaptive_'+"_".join(map(str,['N', N,'L',L, 'Nsigma', Nsigma,  'Nr', Nr]))+'*pickle'
file_list = glob.glob(fmask)

if (len(file_list)==0):
    print "no such files", fmask
    exit()

for fi,fname in enumerate(file_list):
	print 'load', Nsigma,Nr,L
	try:
		with open(fname, 'rU') as f:
			data = pickle.load(f)

		T2.extend(data['T2'])
		T3.extend(data['T3'])
		T4.extend(data['T4'])
		TMRCA.extend(data['TMRCA'])
		diversity.extend(data['selected_diversity'])
		pair_freq.extend(data['LD'])
		
		if fi>0:
			SFS[2]+=np.ma.masked_invalid(data['SFS'][2]).filled(fill_value=0)
			SFS[3]+=np.ma.masked_invalid(data['SFS'][3]).filled(fill_value=0)
		else:
			SFS=list(data['SFS'])
			SFS[2]=np.ma.masked_invalid(SFS[2]).filled(fill_value=0)
			SFS[3]=np.ma.masked_invalid(SFS[3]).filled(fill_value=0)
	except:
		print "Error loading", fname
		pass

pair_freq=np.asarray(pair_freq)

nbins = 15
LD = np.zeros(nbins)
LD_bins = np.zeros(nbins)


print 'process', label

SFS[2]/=len(T2)
SFS[3]/=len(T2)
tmproots = np.roots([1, label[2], -label[1]**2])
sigma = np.real(np.min(tmproots[np.where(np.isreal(tmproots)) and np.where(tmproots>0)]))
length = 1.0/(1.0+label[2]*np.mean(T2)/label[0])

bins = np.arange(0,nbins+1)*int(label[-1]*length/5)
LD_bins=0.5*(bins[1:]+bins[:-1])
for bi in range(nbins):
    ii = np.where((pair_freq[:,0]<bins[bi+1])*(pair_freq[:,0]>bins[bi])*(pair_freq[:,2]+pair_freq[:,3]>0.4))[0]
    LD[bi]=np.mean(LD_func(pair_freq[ii,1], pair_freq[ii,2], pair_freq[ii,3]))


with open('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_new/processed_data_'
          +"_".join(map(str,['N', N,'L',L, 'Nsigma', Nsigma, 'Nr', Nr]))+'.pickle', 
          'w') as f:
    tmp_data = {
        'label':list(label),
        'T2':[np.mean(T2), np.var(T2)],
        'T3':[np.mean(T3), np.var(T3)],
        'T4':[np.mean(T4), np.var(T4)], 
        'TMRCA': [np.mean(TMRCA), np.var(TMRCA)], 
        'sigma':sigma, 'length':length, 'LD':[LD_bins, LD],
        'SFS':SFS}
        
    pickle.dump(tmp_data, f)

