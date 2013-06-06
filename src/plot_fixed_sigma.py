import numpy as np
from matplotlib import pyplot as plt
from Bio import Phylo as P
from matplotlib import cm
import pickle
import glob
import BSC
params = {'backend': 'ps',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
plt.rcParams.update(params)

file_list=glob.glob('../data_fixed_sigma/processed_data*.pickle')
params=[]
mT2=[]
mT3=[]
mT4=[]
mTMRCA=[]
sigma=[]
length=[]
LD_bins=[]
LD=[]
SFS=[]
ni=0
si=2
ri=3
li=1
blocks = []

for fname in file_list:
    with open(fname, 'r') as f:
        tmp_data = pickle.load(f)
    params.append(tmp_data['params'])
    mT2.append(tmp_data['T2'])
    mT3.append(tmp_data['T3'])
    mT4.append(tmp_data['T4'])
    mTMRCA.append(tmp_data['TMRCA'])
    sigma.append(tmp_data['sigma'])
    length.append(tmp_data['length'])
    LD_bins.append(tmp_data['LD'][0])
    LD.append(tmp_data['LD'][1])
    SFS.append(tmp_data['SFS'])
    blocks.append([tmp_data['block_bins'], tmp_data['block_dis']])

params=np.asarray(params)
mT2=np.asarray(mT2)
mT3=np.asarray(mT3)
mT4=np.asarray(mT4)
mTMRCA=np.asarray(mTMRCA)
Nsigma=params[:,ni]*np.real(np.asarray(sigma))
length=np.real(np.asarray(length))
LD_bins=np.asarray(LD_bins)
LD=np.asarray(LD)

T2_length = params[:,li]/(1+params[:,ri]*mT2[:,0])

div=np.zeros_like(mT2[:,0])
for i in range(len(mT2)):
    div[i]=(np.sum(SFS[i][0]*(1-SFS[i][0])*SFS[i][3])/params[i,li])


#ii = np.arange(len(mT2)) 
LD_ii = np.where((T2_length<0.3*params[:,li])*(T2_length>30))[0]
ii=np.where((T2_length<0.3*params[:,li]))[0]
#ii = np.where((div*T2_length>1.5)*(1.0*T2_length/params[:,li]<0.3))[0]


######################
# Figure T2
######################
plt.figure()
plt.scatter((Nsigma[ii]/np.sqrt(2*np.log(Nsigma[ii]+1))), (mT2[ii,-2]/params[ii,ni]), c=np.log10(T2_length[ii]/params[ii,li]), s=50)
x=np.logspace(-1.5, 1.5, 20)
plt.plot([0.1, 10], [1,1], lw=2, ls='-', c='k', label=r'$1$')
plt.plot([1,400], [4, 0.01], lw=2, ls='--', c='k', label=r'$c\sqrt{2\log N\sigma_b}/(N\sigma_b)$')
plt.text(4000, 7, r'$\xi_b/L$', fontsize=24)
logticks = range(-3,1)
cbar =plt.colorbar(ticks=logticks)
cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in logticks])
plt.xlabel(r'$N\sigma_b/\sqrt{2\log N\sigma_b}$')
plt.ylabel(r'$T_2/N$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(loc=3)
plt.savefig('../figures/T2.pdf', bbox_inches='tight', pad_inches=0.4)
plt.savefig('../figures/T2.svg')




######################
# Figure TMRCA
######################
plt.figure()
plt.scatter((Nsigma[ii]/np.sqrt(2*np.log(1+Nsigma[ii]))), (mTMRCA[ii,-2]/params[ii,ni]), c=np.log10(params[ii,ri]+0.01), s=50)
plt.plot(x, 2.0/(1+x))
plt.colorbar()
plt.xlabel(r'$N\sigma/\sqrt{2\log N\sigma}$')
plt.ylabel(r'$T_{MRCA}/N$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.savefig('../figures/TMRCA.pdf')
plt.savefig('../figures/TMRCA.svg')


######################
# Figure T3
######################
# plt.figure()
# plt.scatter((Nsigma[ii]/np.sqrt(2*np.log(1+Nsigma[ii]))), (mT3[ii,-2]/params[ii,ni]), c=np.log10(params[ii,ri]+0.01), s=50)
# x=np.logspace(-1.5, 1.5, 20)
# plt.plot(x, 4.0/3.0/(1+x))
# plt.colorbar()
# plt.xlabel(r'$N\sigma/\sqrt{2\log N\sigma}$')
# plt.ylabel(r'$T_3/N$')
# ax=plt.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.savefig('../figures/T3.pdf')

# ######################
# # Figure T4
# ######################
# plt.figure()
# plt.scatter((Nsigma[ii]/np.sqrt(2*np.log(1+Nsigma[ii]))), (mT4[ii,-2]/params[ii,ni]), c=np.log10(params[ii,ri]+0.01), s=50)
# x=np.logspace(-1.5, 1.5, 20)
# plt.plot(x, 3.0/2.0/(1+x))
# plt.colorbar()
# plt.xlabel(r'$N\sigma/\sqrt{2\log N\sigma}$')
# plt.ylabel(r'$T_4/N$')
# ax=plt.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.savefig('../figures/T4.pdf')

####################
# length vs parameters
###################
# plt.figure()
# plt.scatter(length[ii],params[ii,si]**1.5*params[ii,ri]**-1.5, 
#             c=np.log10(params[ii,ri]+0.01), s=50)
# plt.colorbar()
# plt.xlabel(r'block length')
# plt.ylabel(r'$\sqrt{\mu s^2}/\rho^{3/2}$')
# ax=plt.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.savefig('../figures/length.pdf')

####
# LD vs Distance
####
cset = np.log10(T2_length)
cmin=np.log10(30)
plt.figure()
for i in LD_ii:
    plt.plot(LD_bins[i]/T2_length[i], LD[i]/LD[i][0],
             c=cm.jet((cset[i]-cmin)
                      /(cset.max()-cmin)))
#fake a colorbar...
plt.text(3.4, 1.15, r'$\xi_b$', fontsize=24)
sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.normalize(vmin=cmin, vmax=(max(cset[ii]))))
sm._A = []
tickspos= [30,100,300,1000]
logticks=np.log10(tickspos)
cbar = plt.colorbar(sm, ticks=logticks)
cbar.set_ticklabels([r'$'+str(i)+'$' for i in tickspos])

plt.ylim([0,1.2])
plt.xlabel(r'$d/\xi_b$')
plt.ylabel(r'normalized $r^2$')
plt.savefig('../figures/LD.pdf', bbox_inches='tight', pad_inches=0.2)
plt.savefig('../figures/LD.svg')

############################
# SFS
##########################
BSCN=10000
plt.figure()
BSC_bins, BSC_SFS = BSC.get_BSC_SFS()
for i in range(len(SFS)):
    plt.plot(np.log(SFS[i][0]/(1-SFS[i][0])), 
             SFS[i][2]/SFS[i][1]/params[i,ni]**1, 
             c=cm.jet((np.log10(Nsigma[i])-np.log10(Nsigma.min()))
                      /(np.log10(Nsigma.max())-np.log10(Nsigma.min()))))

#fake a colorbar...
plt.text(10, 5e4, r'$N\sigma_b$', fontsize=24)
sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.normalize(vmin=np.log10(min(Nsigma[ii])), vmax=np.log10(max(Nsigma[ii]))))
sm._A = []
logticks = range(-1,4)
cbar = plt.colorbar(sm, ticks=logticks)
cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in logticks])

plt.xlabel(r'$\rm{derived\ allele\ frequency}\ \nu$')
nutick = np.array([1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1-1e-1, 1-1e-2, 1-1e-3, 1-1e-4])
nutick_labels = []
nutick_labels.extend([r'$0.'+"".join(['0']*(i-1))+'1$' for i in range(4,0,-1)])
nutick_labels.append(r'$0.5$')
nutick_labels.extend([r'$0.'+"".join(['9']*i)+'$' for i in range(1,5)])
plt.xticks(np.log(nutick[1::2]/(1-nutick[1::2])), nutick_labels[1::2])
plt.plot(np.log(SFS[0][0]/(1-SFS[0][0])), 2.0/SFS[0][0], c='k', lw=3, label='Kingman')
plt.plot(np.log(BSC_bins[BSCN]/(1-BSC_bins[BSCN])), BSC_SFS[BSCN]/BSC_SFS[BSCN][-6],  c='k', lw=3,ls='--', label='BSC')
ax=plt.gca()
ax.set_yscale('log')
plt.xlim([-8,8])
plt.legend(loc=1)
plt.savefig('../figures/SFS.pdf', bbox_inches='tight', pad_inches=0.2)         
plt.savefig('../figures/SFS.svg')         


############################
# SFS
##########################
plt.figure()
plt.title(r'color = $\log_{10}(N\sigma_b)$')
for i in range(len(SFS)):
    plt.plot(np.log(SFS[i][0]/(1-SFS[i][0])), 
             SFS[i][3]/SFS[i][1]/params[i,li], 
             c=cm.jet((np.log10(Nsigma[i])-np.log10(Nsigma.min()))
                      /(np.log10(Nsigma.max())-np.log10(Nsigma.min()))))
#fake a colorbar...
sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.normalize(vmin=np.log10(min(Nsigma[ii])), vmax=np.log10(max(Nsigma[ii]))))
sm._A = []
plt.colorbar(sm)

plt.xlabel(r'$\rm{derived\ allele\ frequency}\ \nu$')
nutick = np.array([1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1-1e-1, 1-1e-2, 1-1e-3, 1-1e-4])
nutick_labels = []
nutick_labels.extend([r'$0.'+"".join(['0']*(i-1))+'1$' for i in range(4,0,-1)])
nutick_labels.append(r'$0.5$')
nutick_labels.extend([r'$0.'+"".join(['9']*i)+'$' for i in range(1,5)])
plt.xticks(np.log(nutick[1::2]/(1-nutick[1::2])), nutick_labels[1::2])
plt.plot(np.log(SFS[0][0]/(1-SFS[0][0])), 2.0/SFS[0][0], c='k', lw=3, label='Kingman')
plt.plot(np.log(BSC_bins[BSCN]/(1-BSC_bins[BSCN])), BSC_SFS[BSCN]/BSC_SFS[BSCN][-6],  c='k', lw=3, label='BSC')
ax=plt.gca()
ax.set_yscale('log')
plt.xlim([-6,6])
ax=plt.gca()
ax.set_yscale('log')
plt.savefig('../figures/SFS_sel.pdf')         

###### hap block distributions
plt.figure()
for bl in blocks:
    print bl[1]
    if (len(bl[1])==len(bl[0])):
        plt.plot(bl[0], bl[1])

