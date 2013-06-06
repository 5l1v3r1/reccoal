import numpy as np
import FFPopSim as h
from Bio import Phylo as P
import pop_pickle
import argparse


#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a purifying selection with recombination")
parser.add_argument('--pop', default=1000, type=int, help='Population size (N)')
parser.add_argument('--L', default=1000, type=int, help='number of loci (L)')
parser.add_argument('--r', default=0.0, type=float, help='crossover rate (Nrho)')
parser.add_argument('--sigma', default=0.05, type=float, help='fitness std dev (N\sigma)')
params=parser.parse_args()
N=params.pop
L=params.L
s=0.001
sigma = params.sigma
rho=params.r/L

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


############
# initial buring before variance is pinned
###########
mytree = pop.genealogy.get_tree(gen_loci[0])
MRCA_time = mytree.MRCA.age
while MRCA_time<max(100, pop.generation-MRCA_time):
    pop.evolve()
    if pop.generation%100==0:
        mytree = pop.genealogy.get_tree(gen_loci[0])
        MRCA_time = mytree.MRCA.age
        print "generation:",pop.generation, "until the MRCA has moved past gen/2. current ", MRCA_time

################
# further burnin with pinned variance
###############
switch_time=pop.generation
while pop.generation-switch_time<10*max(100,pop.generation-MRCA_time):
    pop.evolve()
    actual_sigma = np.sqrt(pop.get_fitness_statistics().variance)
    pop.trait_weights = pop.trait_weights*sigma/actual_sigma
    if pop.generation%100==0:
        mytree = pop.genealogy.get_tree(gen_loci[0])
        MRCA_time = mytree.MRCA.age
        print "generations since switch:",pop.generation-switch_time, "until the MRCA has moved past 10 times the typical TMRCA. current TMRCA sits at ", pop.generation - MRCA_time

fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/seeds_fixed_sigma/pop_'+"_".join(map(str,['N', N,'L',L, 'sigma', np.round(sigma,4),  'r', np.round(rho*L, 4)]))+'.pickle'
pop_pickle.dump_pop(pop, fname)
