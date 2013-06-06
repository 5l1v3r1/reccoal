import os
import numpy as np
import glob

def parse_file_name(fname):
    entries = fname.split('_')
    params={}
    for i in range(len(entries)-1, 0, -1):
        try:
            params[entries[i-1]]=float(entries[i])
        except:
            pass
    return params

seed_dir = 'seeds_fixed_sigma'
script = 'run_fixed_sigma.py'
#seed_dir = 'seeds_fixed_effect'
#script = 'run_fixed_effect.py'
#seed_dir = 'seeds_dynbalance'
#script = 'run_dynbalance.py'

nsamples = 30

mem = {300:'1G', 1000:'1G', 3000:'3G', 10000:'8G'}
rt = {300:'0:59:0', 1000:'0:59:0', 3000:'5:59:0', 10000:'9:59:0'}

seed_list = glob.glob('/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/'+seed_dir+'/pop*pickle')
print len(seed_list)
rid=10
for rid in range(11, 20):
    for fname in seed_list:
        params = parse_file_name(fname[:-7])
        params['N']=int(params['N'])    
        params['L']=int(params['L'])
        arg_list = ['--pop', params['N'],'--r', params['r'], '--L', params['L']]
        if 'sigma' in params: arg_list.extend(['--sigma', params['sigma']])
        elif 's' in params: arg_list.extend(['--s', params['s']])
        if 'U' in params: arg_list.extend(['--U', params['U']])
        arg_list.extend(['--n', nsamples, '--id', rid])
        arguments = " ".join(map(str, arg_list))
        #print params
        if (params['N']==3000 and params['L']==3000 and 
            params['sigma']==0.03 and (params['r']>0.3 or params['r']<0.03)):
            command='qsub -cwd -l h_vmem='+mem[params['N']]+' -l h_rt='+rt[params['N']]+'  submit_script.py ' +script+' '+arguments
            os.system(command)
            print command

