__author__ = 'Sancar'

from gurobipy import *
import re
import os
import time
import copy

import numpy as np
import ggplot as gg
from scipy import stats as st
from pandas import DataFrame as pd_df


os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
import pymatlab



pr_dir = os.environ.get('PROJECT_DIR')
session = pymatlab.session_factory(options=" -logfile sgm_ip_run_log"+time.strftime("%H:%M:%S")+".txt")

#pr_dir = 'G:\Sancar\Documents\projects'
pr_dir = '~/projects/'
print pr_dir

pr_dir =  '/home/sadali1/projects/'
root_dir  = os.path.join(pr_dir, 'SeededGraphMatch')
models_dir = os.path.join(root_dir, 'optim_models')


chdir_script = """cd """ + root_dir
#chdir_script = """cd G:\\Sancar\\Documents\\projects\\SeededGraphMatch"""
session.putvalue('chdir_script', chdir_script)


session.run('eval(chdir_script)')
session.run('cur_dir  = pwd')
returndir = session.getvalue('cur_dir')
print returndir
session.run('addpath(genpath(pwd))')
session.run("addpath('lib')")
session.run("addpath('src')")


truematch_orig = []
runtime_orig = []
truematch_mod = []
runtime_mod = []

subg_sizes = [118, 83, 78]

nmc = 0
subgraph_id  = []

for fname in os.listdir(models_dir):
    if fname.endswith(".lp"):
        try:
            pathfile = os.path.join(models_dir, fname)
            print "Processing" + fname
            gmodel = read(pathfile)
            gmodel.update()
            gmodel_relaxed = gmodel.copy()
            allConstr = gmodel_relaxed.getConstrs()
            fname_ints = re.findall(r'\d+', fname)
            subg_index = int(fname_ints[0])
            subgraph_id.append(subg_index)
            mplusn = subg_sizes[subg_index-1];
            m = int(fname_ints[1])
            n= mplusn -m
            mc = int(fname_ints[2])
            relaxedConstr_22= allConstr[0:(n*n)]
            relaxedConstr_12= allConstr[(n*n):(n*n+2*m*n)]
            pen_Constr_22 = [1.0]*len(relaxedConstr_22)
            pen_Constr_12 = [5.0]*len(relaxedConstr_12)
            gmodel_relaxed.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, pen_Constr_22+pen_Constr_12)
            gmodel_relaxed.update()
            gmodel_relaxed.write(pathfile+'relaxed.lp')
            optresults = gmodel_relaxed.optimize()
            print("Relaxed problems solved")
            v_alt = gmodel_relaxed.getVars()
            P_alt = np.zeros((n*n, 1))
            testm = n
            testm_npy = np.array([testm])
            session.putvalue('testm', testm_npy)
            mplusn_npy = np.array([mplusn])
            session.putvalue('mplusn', mplusn_npy)
            seed_m_npy = np.array([m])
            session.putvalue('seedm', seed_m_npy)
            print ("MATLAB vars created")
            session.run("true_perm = 1:mplusn")
            session.run("true_perm= double(true_perm)")
            session.run("seedm= double(seedm)")
            print ("MATLAB vars converted to double")
            it = 0
            for var_a in v_alt[0:(n*n)]:
                P_alt[it, 0] = var_a.X
                it += 1
            P_alt = P_alt.reshape((n, n)).transpose()
            session.putvalue('mod_Phat', P_alt)
            run_script = '[mod_tm_m'+str(mc) +', mod_delta_m'+str(mc) +'] = PermMat2TrueMatch(mod_Phat, true_perm, seedm)'
            try:
                session.run(run_script)
                truematch_mod.append(session.getvalue("mod_delta_m"+str(mc)))
                runtime_mod.append(gmodel_relaxed.Runtime)
            except RuntimeError:
                print ("unable to evaluate true matchin perf of gurobi model in "+fname)
        except (OSError, NameError, ValueError,RuntimeError):
            print "unable to process gurobi model in "+fname



np.savetxt("runtime.csv", runtime_mod, delimiter=",")
np.savetxt("tm.csv", runtime_mod, delimiter=",")

tm_df = pd_df({
    "tm": truematch_orig + truematch_mod,
    "runtime": runtime_orig + runtime_mod,
    "orig_or_mod": ["mod"]*len(runtime_mod) # ["orig"]*nmc
})
nmc = len(runtime_mod)




#print gg.ggplot(tm_df, aes('orig_or_mod', 'runtime')) + \
#  gg.geom_line(colour='steelblue')

comp_plot = gg.ggplot(data=tm_df, aesthetics=gg.aes(x='runtime', y='tm')) + gg.geom_point() + gg.scale_x_log10()
gg.ggsave("graphmatch_IP_runtime_vs_tm.pdf",plot = comp_plot)
print   "mean of tm:"+ str(np.mean( truematch_mod))
