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
s=params.s
r=params.r
U=params.U

fixed_sigma=0
fixed_effect=1
dynbalance=2

T2=[]
T3=[]
T4=[]
pair_blocks = []
TMRCA=[]
diversity=[]
fitness_variance = []
pair_freq=[]

arg_list = ['N',N, 'L', L]
parameter_list = [N,L]
if sigma>0:
    simtype=fixed_sigma
    data_dir='data_fixed_sigma'
    arg_list.extend(['sigma', sigma])
    parameter_list.append(sigma)
elif s>0 and U<0:
    simtype=fixed_effect
    data_dir='data_fixed_effect'
    arg_list.extend(['s', s])
    parameter_list.append(s)
elif U>0 and sigma<0:
    simtype=dynbalance
    data_dir='data_dynbalance'
    arg_list.extend(['s', s, 'U', U])
    parameter_list.extend([s,U])
arg_list.extend(['r', r])
parameter_list.append(r)

fmask = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/'+data_dir+\
    '/hist_adaptive_'+"_".join(map(str,arg_list))+'*pickle'
file_list = glob.glob(fmask)

if (len(file_list)==0):
    print "no such files", fmask
    exit()

for fi,fname in enumerate(file_list):
    print 'load', fname
    try:
        with open(fname, 'rU') as f:
            data = pickle.load(f)

        T2.extend(data['T2'])
        T3.extend(data['T3'])
        T4.extend(data['T4'])
        try:
            fitness_variance.extend(data['sigmasq'])
        except:
            print "no total variance"
            fitness_variance.append(0)
            pass
        try:
            pair_blocks.extend(data['blocks'])
            print len(pair_block)
        except:
            print "no block length"
            #pair_blocks.append(0)
            pass
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

SFS[2]/=len(T2)
SFS[3]/=len(T2)
if simtype==fixed_sigma:
    tmproots = np.roots([1, r, -sigma**2])
elif simtype==fixed_effect:
    tmproots = np.roots([1, r, -np.mean(fitness_variance)])
elif simtype==dynbalance:
    tmproots = np.roots([1, r,0, -U*s**2])
else:
    adfsdf
    
sigma_sol = np.real(np.min(tmproots[np.where(np.isreal(tmproots)) and np.where(tmproots>0)]))
length = 1.0/(1.0+r*np.mean(T2))

nbins = 15
LD = np.zeros(nbins)
bins = np.arange(0,nbins+1)*int(L*length/5)
LD_bins=0.5*(bins[1:]+bins[:-1])
for bi in range(nbins):
    ii = np.where((pair_freq[:,0]<bins[bi+1])*(pair_freq[:,0]>bins[bi]))[0] #*(pair_freq[:,2]+pair_freq[:,3]>0.4))[0]
    LD[bi]=np.mean(LD_func(pair_freq[ii,1], pair_freq[ii,2], pair_freq[ii,3]))

nbins=50
bins = np.arange(0,nbins+1)*int(L*length/10)
block_bins=0.5*(bins[1:]+bins[:-1])
if (len(pair_blocks)>0):
    block_dis = np.histogram(pair_blocks, bins=bins)[0]
else:
    block_dis=[]


with open('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/'+data_dir+'/processed_data_'
          +"_".join(map(str,arg_list))+'.pickle', 
          'w') as f:
    tmp_data = {
        'params':list(parameter_list),
        'T2':[np.mean(T2), np.var(T2)],
        'T3':[np.mean(T3), np.var(T3)],
        'T4':[np.mean(T4), np.var(T4)], 
        'sigmasq_tot':np.mean(fitness_variance),
        'TMRCA': [np.mean(TMRCA), np.var(TMRCA)], 
        'sigma':sigma_sol, 'length':length, 'LD':[LD_bins, LD],
        'SFS':SFS, 'block_dis':block_dis, 'block_bins':block_bins}
        
    pickle.dump(tmp_data, f)

