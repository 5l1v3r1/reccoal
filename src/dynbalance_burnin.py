import numpy as np
import FFPopSim as h
from Bio import Phylo as P
import pop_pickle
import argparse


#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a purifying selection with recombination")
parser.add_argument('--pop', default=1000, type=int, help='Population size (N)')
parser.add_argument('--L', default=1000, type=int, help='number of loci (L)')
parser.add_argument('--r', default=0.0, type=float, help='crossover rate')
parser.add_argument('--s', default=0.05, type=float, help='fitness effect of individual loci')
parser.add_argument('--U', default=0.05, type=float, help='mutation rate, entire chromosome')
params=parser.parse_args()
N=params.pop
L=params.L
s=params.s
rho=params.r/L
mu =params.U/L

pop = h.haploid_highd(L, all_polymorphic=False)
pop.outcrossing_rate=1.0
pop.crossover_rate=rho
pop.carrying_capacity=N
pop.mutation_rate=mu
gen_loci = [L/2]
pop.track_locus_genealogy(gen_loci)
sel_coeff = -np.ones(L)*s*0.5
sel_coeff[gen_loci]=1e-10
pop.set_fitness_additive(sel_coeff)
pop.set_wildtype(N)



mytree = pop.genealogy.get_tree(gen_loci[0])
MRCA_time = mytree.MRCA.age

################
#  burnin with fixed effect and fluctuating variance
###############
switch_time=pop.generation
fit_stat = []
while pop.generation-switch_time<10*max(100,pop.generation-MRCA_time):
	pop.evolve()
	if pop.generation%100==0:
		mytree = pop.genealogy.get_tree(gen_loci[0])
		MRCA_time = mytree.MRCA.age
		fs = pop.get_fitness_statistics()
		fit_stat.append([pop.generation, fs.mean, fs.variance])
		print "generations since switch:",pop.generation-switch_time, "until the MRCA has moved past 10 times the typical TMRCA. current TMRCA sits at ", pop.generation - MRCA_time

fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/seeds_dynbalance/pop_'+"_".join(map(str,['N', N,'L',L, 's', np.round(s,4), 'U', np.round(mu*L,4), 'r', np.round(rho*L, 4)]))+'.pickle'
pop_pickle.dump_pop(pop, fname)
fitstat_fname = '/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/seeds_dynbalance/fitstat_'+"_".join(map(str,['N', N,'L',L, 's', np.round(s,4), 'U', np.round(mu*L,4), 'r', np.round(rho*L, 4)]))+'.pickle'
np.savetxt(fitstat_fname, fit_stat)
