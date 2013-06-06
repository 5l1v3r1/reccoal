import numpy as np
import FFPopSim as h
from matplotlib import pyplot as plt
from treehelpers import *
from Bio import Phylo as P
import pickle
import argparse
import random as rd
import pop_pickle

#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a purifying selection with recombination")
parser.add_argument('--pop', default=1000, type=int, help='Population size (N)')
parser.add_argument('--L', default=1000, type=int, help='number of loci (L)')
parser.add_argument('--r', default=0.0, type=float, help='crossover rate rho')
parser.add_argument('--s', default=0.05, type=float, help='fitness effect ')
parser.add_argument('--U', default=1e-1, type=float, help='genome wide mutation rate')
parser.add_argument('--n', default=200, type=int, help='number of samples')
parser.add_argument('--id', default='', type=str, help='job_id')

params=parser.parse_args()
N=params.pop
L=params.L
s = params.s
rho=params.r/L
U=params.U
mu=U/L

#time between sampling intervals, specific for dynbalance
dt=int(min(0.5*N,np.sqrt(L*rho/(U*s**2))+100))
print "sampling interval",dt

nsamples = params.n
if (params.id!=''): id='_'+params.id
else: id=''

gen_loci = [L/2]
#load population from seed.
fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/seeds_dynbalance/pop_'+"_".join(map(str,['N', N,'L',L, 's', np.round(s,4), 'U', np.round(U,4),  'r', np.round(rho*L, 4)]))+'.pickle'
pop=pop_pickle.load_pop(fname, gen_loci)

popsample = 100
nT2=nT3=nT4 = 5
nLD =100

#set up bins non-uniformly make histograms of the allele frequencies
sfsbins = 23
bins=np.zeros(sfsbins)
bins[1:-1]=np.exp(np.linspace(-np.log(1.2*N),np.log(1.2*N),sfsbins-2)) \
            /(1+np.exp(np.linspace(-np.log(1.2*N),np.log(1.2*N),sfsbins-2)))
bins[0]=0
bins[-1]=1
bincenters = 0.5*(bins[1:]+bins[:-1])
dx= bins[1:]-bins[:-1]
neutralSFS = np.zeros_like(bincenters)
selectedSFS = np.zeros_like(bincenters)

ancestor_fitness = []
pop_fitness = []
mean_fitness = []
fitness_variance = []
diversity = []
TMRCA = []
pair_freqs = []
T2=[]
T3=[]
T4=[]

mfit_scratch = np.zeros(10000)

###########################
# function to save the collected data to file. uses lots of globals
##########################
def save_data():
    print "sample",n,"out of", nsamples
    diversity_tmp = np.asarray(diversity)
    fitnessbins = np.linspace(-30*s, 30*s,201)

    #produce histograms
    popfithist, x = np.histogram(pop_fitness, fitnessbins, normed=True)
    ancestorfithist, x = np.histogram(ancestor_fitness, fitnessbins, normed=True)

    #dump the results to file
    data_fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/data_dynbalance/hist_adaptive_'\
    +"_".join(map(str,['N', N,'L',L, 's', np.round(s,4), 'U', np.round(U,4),  'r', np.round(rho*L, 4)]))+id+ '.pickle'
    print "saving to", data_fname
    with open(data_fname, 'w') as f:
        pickle.dump({'fithist': (fitnessbins, popfithist, ancestorfithist), 
                     'SFS': (bincenters, dx, neutralSFS, selectedSFS), 
                     'selected_diversity': diversity_tmp, 'sigmasq':np.asarray(fitness_variance),
                     'TMRCA':np.asarray(TMRCA), 'T2':np.asarray(T2), 
                     'T3':np.asarray(T3) ,'T4':np.asarray(T4), 
                     'LD': np.array(pair_freqs)}, f)


################
# further burnin with pinned variance
###############
mytree = pop.genealogy.get_tree(gen_loci[0])
MRCA_time = mytree.MRCA.age
while pop.generation < 2*max(100,pop.generation-MRCA_time):
    pop.evolve()
    if pop.generation%100==0:
        mytree = pop.genealogy.get_tree(gen_loci[0])
        MRCA_time = mytree.MRCA.age
        print "generation:",pop.generation, "until the TMRCA is less than half the generation. current ", pop.generation-MRCA_time

    mfit_scratch[pop.generation%mfit_scratch.shape[0]]=pop.get_fitness_statistics().mean

##############
# measure
##############
for n in xrange(nsamples):
    if (n%3==2): 
        print pop.generation
        save_data()
        
    #advance to next sample
    for ti in xrange(dt):
        pop.evolve()
        mfit_scratch[pop.generation%mfit_scratch.shape[0]]=pop.get_fitness_statistics().mean
    

    # LD and pair frequencies
    daf = pop.get_derived_allele_frequencies()
    af = pop.get_allele_frequencies()
    ii = np.where((af>0.1)*(af<0.9))[0]    
    for i in xrange(min(nLD, len(ii)*(len(ii)-1)/2)):
        locus1=ii[rd.randrange(0,len(ii))]
        locus2=ii[rd.randrange(0,len(ii))]
        pf = pop.get_pair_frequency(locus1,locus2)
        pair_freqs.append([np.abs(locus1-locus2), pf, af[locus1], af[locus2]])

    y,x = np.histogram(daf,  bins=bins)
    selectedSFS+=y
   
    # fitness distribution
    sample = pop.random_clones(popsample)
    mfit= mfit_scratch[pop.generation%mfit_scratch.shape[0]]
    pop_fitness.extend([pop.get_fitness(i)-mfit for i in sample])
    fit_stat = pop.get_fitness_statistics()
    fitness_variance.append(fit_stat.variance)

    # genealogies, coalescent time, and SFS
    for locus in gen_loci:
        T=pop.genealogy.get_tree(locus)
        A=T.nodes
        B=A.keys()
        ancestor_fitness.append(A[B[B.index(T.MRCA)]].fitness-mfit_scratch[T.MRCA.age%mfit_scratch.shape[0]])
        tmpn, tmpSFS = get_SFS(T.to_Biopython_tree())
        y,x = np.histogram(tmpSFS[:,0], weights = tmpSFS[:,1], bins=bins)
        neutralSFS+=y
        TMRCA.append(pop.generation-T.MRCA.age)
 
        tmpT2 = []
        for i in xrange(nT2):
            sam2 = rd.sample(T.leafs,2)
            tmpT = T.create_subtree_from_keys(sam2)
            tmpT2.append(pop.generation-tmpT.MRCA.age)
        T2.append(tmpT2)

        tmpT3 = []
        if (len(T.leafs)>3):
            for i in xrange(nT3):
                sam3 = rd.sample(T.leafs,3)
                tmpT = T.create_subtree_from_keys(sam3)
                tmpT3.append(pop.generation-tmpT.MRCA.age)
            T3.append(tmpT3)

        tmpT4 = []
        if (len(T.leafs)>4):
            for i in xrange(nT4):
                sam4 = rd.sample(T.leafs,4)
                tmpT = T.create_subtree_from_keys(sam4)
                tmpT4.append(pop.generation-tmpT.MRCA.age)
            T4.append(tmpT4)

    tmpdiv = pop.get_diversity_statistics(100)
    diversity.append([tmpdiv.mean, tmpdiv.variance])

#############
# final saving of data
##############
save_data()

#################
# replace seed population
################
fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/seeds_dynbalance/pop_'\
    +"_".join(map(str,['N', N,'L',L, 's', np.round(s,4),  'U', np.round(U,4), 'r', np.round(rho*L, 4)]))+'.pickle'
pop_pickle.dump_pop(pop, fname)


