import numpy as np
import pickle
import glob

def get_BSC_SFS():
	BS_af = {}
	BS_af_bins = {}
	BS_file_count = {}
	'''LOAD SFS OF DIRECT BSC SIMULATIONS'''
	filelist = glob.glob('/ebio/ag-neher/share/users/rneher/Coalescent/Simulations/BS_direct2/SFS*.pickle')
	for fname in filelist:
		entries=fname.split('_')
		N=int(entries[-2])
		file=open(fname, 'r')
		temp, temp_bins=pickle.load(file)
		if (N in BS_af):
			BS_af[N]+=temp/(BS_af_bins[N][1:]-BS_af_bins[N][:-1])
			BS_file_count[N]+=1
		else:
			BS_af_bins[N]=temp_bins
			BS_file_count[N]=1
			BS_af[N]=temp/(BS_af_bins[N][1:]-BS_af_bins[N][:-1])
	#Produce centered bins
	for N in BS_af_bins:
		BS_af_bins[N]=0.5*(BS_af_bins[N][1:]+BS_af_bins[N][:-1])
	#divide the spectra by the file count to normalize
	for N in BS_af:
		BS_af[N]/=BS_file_count[N]

	return BS_af_bins, BS_af
