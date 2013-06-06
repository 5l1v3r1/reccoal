import numpy as np

def get_mut_count_subtree(R):
	'''recursive function that returns a list of pairs of branch length and the number of down stream leaves\
	branch length is proportional to the number of mutations, the number of leaves to the frequency of the mutations.'''
	if R.is_terminal():
		return int(R.name.split('_')[1]),[]				#a terminal branch has no downstream mutations
	else:
		tmp_mut_count = []
		sub_tree_size = 0
		for C in R.clades:		
			#loop over all children of the node and accumulate their mutations
			csize,cmut  = get_mut_count_subtree(C)
			sub_tree_size += csize
			tmp_mut_count.extend(cmut)
			#add mutations that happened on the branch to the child. those are present in csize leafs
			tmp_mut_count.append([csize,C.branch_length])
		return sub_tree_size, tmp_mut_count


def get_SFS(T):
	'''returns the sample size and the site frequency spectrum of a tree T'''
	#get a list of all opportunities for mutations on the tree below the root.
	#all derived mutations happen there
	sample_size,SFS = get_mut_count_subtree(T.root)
	#convert to a numpy array and normalize
	SFS = np.asarray(SFS)
	SFS[:,0]/=sample_size
	return sample_size,SFS
