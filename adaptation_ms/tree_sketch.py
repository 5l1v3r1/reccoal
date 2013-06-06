import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd
import Bio
from Bio import Phylo

print "This script is meant to produce to example plots of trees, with and without selection"


L = 1000    #number of segregating sites
s = 1e-2    #single site effect
N = 10000   #population size
r = 0.0     #outcrossing rate
sigma = 0.05
sigmasq=sigma**2

sample_size=30  #number of individuals whose genealogy is looked at
nsamples = 3
burnin = 2000   #either ~5*N or 5/s, depending on whether coalescence is dominated by drift or draft
dt = 1000       #time between samples

#set up population, switch on infinite sites mode
pop=h.haploid_highd(L,all_polymorphic=True)

#set the population size via the carrying capacity
pop.carrying_capacity= N

#set the outcrossing_rate to 0
pop.outcrossing_rate = 0.0

#set the effect sizes of the mutations that are injected (the same at each site in this case)
pop.set_fitness_additive(np.ones(L)*s*0.5)

#track the genealogy at a central locus L/2 (which one doesn't matter in the asexual case)
pop.track_locus_genealogy([L/2])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

print "Population parameters:"
pop.status()

#burn in
print "\nEquilibrate:"
while pop.generation<burnin:
    if (pop.generation%100==0): print "Burn in: at", pop.generation, "out of", burnin, "generations"
    pop.evolve()

while pop.generation<2*burnin:
    if (pop.generation%100==0): print "Burn in: at", pop.generation, "out of", burnin, "generations"
    pop.evolve()
    fit_stat = pop.get_fitness_statistics()
    actual_variance =fit_stat.variance
    pop.trait_weights = np.sqrt(sigmasq/actual_variance)*pop.trait_weights



print "\nPlot coalescent trees:"
fig=plt.figure(figsize=(7,10))
fig.suptitle("".join(map(str,['N=',N,'  r=',r,'  L=',L, '  s=',s])), fontsize=18)
for si in xrange(nsamples):
    print "sample",si,"out of",nsamples
    #evolve a while before sampling the next tree
    for ii in range(dt):
        pop.evolve()
        fit_stat = pop.get_fitness_statistics()
        actual_variance =fit_stat.variance
        pop.trait_weights = np.sqrt(sigmasq/actual_variance)*pop.trait_weights
   
    #draw a sample from the population, convert its genealogy to a BioPython tree object and plot
    tree = pop.genealogy.get_tree(L/2)
    subtree = tree.create_subtree_from_keys(rd.sample(tree.leafs,sample_size)).to_Biopython_tree()
    subtree.ladderize()
    ax = plt.subplot(3,1,si+1)
    if (float(Bio.__version__) >=1.60): 
        Phylo.draw(subtree,label_func=lambda x:"", axes=ax)
        plt.xlim(subtree.root.branch_length-50, pop.generation+50)
    else: Phylo.draw(subtree,label_func=lambda x:"")
    plt.draw()

plt.savefig('../figures/trees_positive_'+"".join(map(str,['N=',N,'_r=',r,'_L=',L, '_sigma=',sigma,'.svg'])))

