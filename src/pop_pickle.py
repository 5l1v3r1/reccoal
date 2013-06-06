import numpy as np
import FFPopSim as h
import pickle
import bz2

def dump_pop(pop, fname):
    '''dumps a population to file, epistasis not supported'''
    pop_dict = {}
    pop_dict['gts'] = pop.get_genotypes()
    pop_dict['N'] = pop.carrying_capacity
    pop_dict['L'] = pop.L
    pop_dict['mu'] = pop.mutation_rate
    pop_dict['xo'] = pop.crossover_rate
    pop_dict['ox'] = pop.outcrossing_rate
    pop_dict['circ'] = pop.circular
    pop_dict['clone_sizes'] = pop.get_clone_sizes()
    pop_dict['recmodel'] = pop.recombination_model
    pop_dict['traits'] = [pop.get_trait_additive(i) for i in range(pop.number_of_traits)]
    pop_dict['all_polymorphic']  = pop.all_polymorphic
    pop_dict['ancestral'] = pop.get_ancestral_states()
    pop_dict['trait_weights'] = pop.trait_weights
    with open(fname, 'wb') as f:
        f.write(pickle.dumps(pop_dict,pickle.HIGHEST_PROTOCOL).encode('bz2'))

def load_pop(fname, gen_loci=[]):
    '''loads a population from a compressed pickle file and restores all parameters'''
    with bz2.BZ2File(fname, 'rb') as f:
        pop_dict = pickle.load(f)

    pop=h.haploid_highd(pop_dict['L'], all_polymorphic = pop_dict['all_polymorphic'], number_of_traits = len(pop_dict['traits']))
    pop.carrying_capacity = pop_dict['N']
    if (pop.all_polymorphic==False): pop.mutation_rate = pop_dict['mu']
    pop.crossover_rate = pop_dict['xo']
    pop.outcrossing_rate = pop_dict['ox']
    pop.circular = pop_dict['circ']

    pop.recombination_model = pop_dict['recmodel']
    for i in range(pop.number_of_traits):
        pop.set_trait_additive(pop_dict['traits'][i], i)
    
    pop.trait_weights=pop_dict['trait_weights']
    if len(gen_loci)>0: pop.track_locus_genealogy(gen_loci)
    pop.set_genotypes_and_ancestral_state(pop_dict['gts'], 
                                          pop_dict['clone_sizes'], 
                                          pop_dict['ancestral'])
    
    return pop







