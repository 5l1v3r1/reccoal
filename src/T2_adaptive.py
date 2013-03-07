import numpy as np
import FFPopSim as h
from matplotlib import pyplot as plt
from treehelpers import *
from Bio import Phylo as P
import pickle
import argparse
import random as rd


#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a purifying selection with recombination")
parser.add_argument('--pop', default=10000, type=int, help='Population size (N)')
parser.add_argument('--L', default=10000, type=int, help='number of loci (L)')
parser.add_argument('--Nrho', default=0.0, type=float, help='Scaled crossover rate (Nrho)')
parser.add_argument('--Nsigma', default=100.0, type=float, help='Scaled fitness std dev (N\sigma)')
parser.add_argument('--dt', default=200, type=int, help='sampling interval')
parser.add_argument('--n', default=200, type=int, help='number of samples')
parser.add_argument('--id', default='', type=str, help='job_id')
params=parser.parse_args()
N=params.pop
L=params.L
s=0.001
sigma = params.Nsigma/N
rho=params.Nrho/N/L
dt=params.dt
nsamples = params.n
if (params.id!=''): id='_'+params.id
else: id=''

pop = h.haploid_highd(L, all_polymorphic=True)
pop.outcrossing_rate=1.0
pop.crossover_rate=rho
pop.carrying_capacity=N
gen_loci = [L/2]
pop.track_locus_genealogy(gen_loci)
sel_coeff = np.ones(L)*s*0.5
sel_coeff[gen_loci]=1e-10
pop.set_fitness_additive(sel_coeff)
pop.set_wildtype(N)

burnin = int(5*N)
popsample = 100
nT2=nT3=nT4 = 5

#set up bins non-uniformly make histograms of the allele frequencies
bins=np.exp(np.linspace(-2.5*np.log(10),2.5*np.log(10),21)) \
			/(1+np.exp(np.linspace(-2.5*np.log(10),2.5*np.log(10),21)))
bins[0]=1.5/pop.carrying_capacity; bins[-1]=1
bincenters = 0.5*(bins[1:]+bins[:-1])
dx= bins[1:]-bins[:-1]
neutralSFS = np.zeros_like(bincenters)

ancestor_fitness = []
pop_fitness = []
mean_fitness = []
diversity = []
TMRCA = []
T2=[]
T3=[]
T4=[]

mfit_scratch = np.zeros(10000)
for ti in xrange(burnin):
    pop.evolve()
    if pop.generation%1000==0:
        print pop.generation, "out of", burnin
    mfit_scratch[pop.generation%mfit_scratch.shape[0]]=pop.get_fitness_statistics().mean

for ti in xrange(burnin):
    pop.evolve()
    actual_sigma = np.sqrt(pop.get_fitness_statistics().variance)
    pop.trait_weights = pop.trait_weights*sigma/actual_sigma
    if pop.generation%1000==0:
        print pop.generation, "out of", burnin
    mfit_scratch[pop.generation%mfit_scratch.shape[0]]=pop.get_fitness_statistics().mean


for n in xrange(nsamples):
    if (n%10==0): print "sample",n,"out of", nsamples
    for ti in xrange(dt):
        pop.evolve()
        actual_sigma = np.sqrt(pop.get_fitness_statistics().variance)
        pop.trait_weights = pop.trait_weights*sigma/actual_sigma
        mfit_scratch[pop.generation%mfit_scratch.shape[0]]=pop.get_fitness_statistics().mean
    mfit= mfit_scratch[pop.generation%mfit_scratch.shape[0]]
    sample = pop.random_clones(popsample)
    pop_fitness.extend([pop.get_fitness(i)-mfit for i in sample])
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

diversity = np.asarray(diversity)
fitnessbins = np.linspace(-10*sigma, 10*sigma,201)
popfithist, x = np.histogram(pop_fitness, fitnessbins, normed=True)
ancestorfithist, x = np.histogram(ancestor_fitness, fitnessbins, normed=True)

with open('/ebio/ag-neher/share/users/rneher/WeakSelectionCoalescent/data/hist_adaptive_'
          +"_".join(map(str,['N', N,'L',L, 'Nsigma', N*sigma,  'Nrho', N*rho*L]))+id+'.pickle', 'w') as f:
    pickle.dump({'fithist': (fitnessbins, popfithist, ancestorfithist), 
                 'SFS': (bincenters, dx, neutralSFS), 
                 'selected_diversity': diversity, 
                 'TMRCA':np.asarray(TMRCA), 'T2':np.asarray(T2), 'T3':np.asarray(T3) ,'T4':np.asarray(T4)}, f)

#fig=plt.figure()
#pdis=plt.hist(pop_fitness, fitnessbins, label='population',normed=True, alpha=.7)
#ancdis=plt.hist(ancestor_fitness, fitnessbins, label='MRCA', normed=True, alpha=.7)
#plt.xlabel('fitness')
#plt.ylabel('distribution')
#plt.xlim([fitnessbins[np.where(pdis[0]>0)[0][0]]+5*s, fitnessbins[np.where(ancdis[0]>0)[0][-1]]-5*s])
#plt.legend(loc=2)
#plt.savefig('dis'+'_'.join(map(str,['N', N, 'Ns', -N*s, 'NU', N*mu*L]))+ '.pdf')

#plt.figure()
#plt.plot(np.log(bincenters/(1.0-bincenters)), neutralSFS/dx)
#ax=plt.gca()
#ax.set_yscale('log')

#P.draw(T.to_Biopython_tree(), label_func = lambda x:"")
#plt.xlim([T.MRCA.age-50, pop.generation+50])
#plt.draw()




