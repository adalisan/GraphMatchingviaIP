#!/usr/bin/python

import sys

__author__ = 'Sancar'

from gurobipy import *
import os
import copy
import pymatlab
import numpy as np
import ggplot as gg
from scipy import stats as st

from pandas import DataFrame as pd_df

pr_dir = os.environ.get('PROJECT_DIR')
session = pymatlab.session_factory()


#pr_dir = '~/projects/'
print pr_dir

root_dir = os.path.join(pr_dir, 'SeededGraphMatch')



n = 35
m = 10

mplusn = m+n
p = 0.2
corr = 0.7
nmc = 20


m_npy = np.array([m])

n_npy = np.array([n])
mplusn_npy = np.array([mplusn])
p_npy = np.array([p])
corr_npy = np.array([corr])



session.putvalue('m', m_npy)
session.putvalue('n', n_npy)
session.putvalue('mplusn', mplusn_npy)
session.putvalue('p', p_npy)
session.putvalue('corr', corr_npy)


convert_to_dbl = [''.join([x, "=double(", x, ")"]) for x in ["m", "n", "mplusn", "p", "corr"]]


for scr in convert_to_dbl:
    print scr
    session.run(scr)

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


for mc in range(nmc):

    run_script = """run ./src/GM_exact_model_general_param.m"""
    #mc_npy = np.array([mc])
    #session.putvalue('tag', mc_npy)
    set_tag =  'tag=uint32('+str(mc)+')'
    print set_tag
    session.run(set_tag)
    print session.getvalue('tag')
    session.putvalue('run_script', run_script)
    session.run('eval(run_script)')
    orig_model_fname = session.getvalue('fname_1')
    infeas_model_fname = session.getvalue('fname_2')
    gmodel = read(os.path.join(returndir, 'src', infeas_model_fname))
    gmodel.update()
    allConstr = gmodel.getConstrs()
    relaxedConstr_22= allConstr[0:(n*n)]
    relaxedConstr_12= allConstr[(n*n):(n*n+2*m*n)]
    pen_Constr_22 = [1.0]*len(relaxedConstr_22)
    pen_Constr_12 = [1.0]*len(relaxedConstr_12)
    gmodel_relaxed = gmodel.copy()
    gmodel_relaxed.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, \
                             pen_Constr_22+pen_Constr_12)
    gmodel_relaxed.update()
    gmodel_relaxed.write(infeas_model_fname + 'relaxed.lp')
    optresults = gmodel_relaxed.optimize()
    v = gmodel_relaxed.getVars()
    P = np.zeros((n*n, 1))
    it = 0
    for var in v[0:(n*n)]:
        #print var.VarName, var.X
        P[it, 0] = var.X
        it += 1
    P = P.reshape((n, n)).transpose()
    session.putvalue('Phat', P)
    #run_script = "Psoln = DS2Perm(Phat)"
    #session.run(run_script)
    run_script = '[tm_m'+str(mc) +', delta_m'+str(mc) +'] = PermMat2TrueMatch(Phat, true_perm, m)'
    session.run(run_script)
    truematch_orig.append(session.getvalue("delta_m"+str(mc)))
    runtime_orig.append(gmodel_relaxed.Runtime)
    pen_Constr_12 = [5.0]*len(relaxedConstr_12)
    gmodel_relaxed_alt = read(os.path.join(returndir, 'src', infeas_model_fname))
    gmodel_relaxed_alt.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, \
                                 pen_Constr_22+pen_Constr_12)
    gmodel_relaxed_alt.update()
    optresults = gmodel_relaxed_alt.optimize()
    v_alt = gmodel_relaxed_alt.getVars()
    P_alt = np.zeros((n*n, 1))
    it = 0
    for var_a in v_alt[0:(n*n)]:
        #print var.VarName, var.X
        P_alt[it, 0] = var_a.X
        it += 1
    P_alt = P_alt.reshape((n, n)).transpose()
    session.putvalue('mod_Phat', P_alt)
    perm_equal = np.allclose(P, P_alt)
    print "Two permutation solutions are same? "+str(perm_equal)
    print "Nonzero elements are "+str(np.nonzero(P-P_alt))
    run_script = '[mod_tm_m'+str(mc) +', mod_delta_m'+str(mc) +'] = PermMat2TrueMatch(mod_Phat, true_perm, m)'
    session.run(run_script)
    truematch_mod.append(session.getvalue("mod_delta_m"+str(mc)))
    runtime_mod.append(gmodel_relaxed_alt.Runtime)
    #tuning = gmodel_relaxed.tune()
    #gmodel_relaxed.getTuneresult(0)
    #gmodel_relaxed.optimize()
    # optimal_runtime[mc] =  gmodel_relaxed.Runtime

tm_df = pd_df({
    "tm": truematch_orig + truematch_mod,
    "runtime": runtime_orig + runtime_mod,
    "orig_or_mod": ["orig"]*nmc + ["mod"]*nmc
})

p_val_truematch_diff = st.ttest_ind(truematch_orig, truematch_mod)
p_val_timediff = st.ttest_ind(runtime_orig, runtime_mod)


#print gg.ggplot(tm_df, aes('orig_or_mod', 'tm')) + \
#  gg.geom_line(colour='steelblue')


#print gg.ggplot(tm_df, aes('orig_or_mod', 'runtime')) + \
#  gg.geom_line(colour='steelblue')







#!/usr/bin/python

import sys

__author__ = 'Sancar'

from gurobipy import *
import os
import copy
import pymatlab
import numpy as np
import ggplot as gg
from scipy import stats as st

from pandas import DataFrame as pd_df

pr_dir = os.environ.get('PROJECT_DIR')
session = pymatlab.session_factory()


#pr_dir = '~/projects/'
print pr_dir

root_dir = os.path.join(pr_dir, 'SeededGraphMatch')



n = 35
m = 10

mplusn = m+n
p = 0.2
corr = 0.7
nmc = 20


m_npy = np.array([m])

n_npy = np.array([n])
mplusn_npy = np.array([mplusn])
p_npy = np.array([p])
corr_npy = np.array([corr])



session.putvalue('m', m_npy)
session.putvalue('n', n_npy)
session.putvalue('mplusn', mplusn_npy)
session.putvalue('p', p_npy)
session.putvalue('corr', corr_npy)


convert_to_dbl = [''.join([x, "=double(", x, ")"]) for x in ["m", "n", "mplusn", "p", "corr"]]


for scr in convert_to_dbl:
    print scr
    session.run(scr)

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

run_script = """run ./src/GM_exact_model_general_param.m"""
for mc in range(nmc):
    #mc_npy = np.array([mc])
    #session.putvalue('tag', mc_npy)
    set_tag ='tag='+str(mc)
    print set_tag
    session.run(set_tag)
    print session.getvalue('tag')
    session.putvalue('run_script', run_script)
    session.run('eval(run_script)')
    orig_model_fname = session.getvalue('fname_1')
    infeas_model_fname = session.getvalue('fname_2')
    gmodel = read(os.path.join(returndir, 'src', infeas_model_fname))
    gmodel.update()
    allConstr = gmodel.getConstrs()
    relaxedConstr_22= allConstr[0:(n*n)]
    relaxedConstr_12= allConstr[(n*n):(n*n+2*m*n)]
    pen_Constr_22 = [1.0]*len(relaxedConstr_22)
    pen_Constr_12 = [1.0]*len(relaxedConstr_12)
    gmodel_relaxed = gmodel.copy()
    gmodel_relaxed.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, \
                             pen_Constr_22+pen_Constr_12)
    gmodel_relaxed.update()
    gmodel_relaxed.write(infeas_model_fname + 'relaxed.lp')
    optresults = gmodel_relaxed.optimize()
    v = gmodel_relaxed.getVars()
    P = np.zeros((n*n, 1))
    it = 0
    for var in v[0:(n*n)]:
        #print var.VarName, var.X
        P[it, 0] = var.X
        it += 1
    P = P.reshape((n, n)).transpose()
    session.putvalue('Phat', P)
    #run_script = "Psoln = DS2Perm(Phat)"
    #session.run(run_script)
    run_script = '[tm_m'+str(mc) +', delta_m'+str(mc) +'] = PermMat2TrueMatch(Phat, true_perm, m)'
    session.run(run_script)
    truematch_orig.append(session.getvalue("delta_m"+str(mc)))
    runtime_orig.append(gmodel_relaxed.Runtime)
    pen_Constr_12 = [5.0]*len(relaxedConstr_12)
    gmodel_relaxed_alt = read(os.path.join(returndir, 'src', infeas_model_fname))
    gmodel_relaxed_alt.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, \
                                 pen_Constr_22+pen_Constr_12)
    gmodel_relaxed_alt.update()
    optresults = gmodel_relaxed_alt.optimize()
    v_alt = gmodel_relaxed_alt.getVars()
    P_alt = np.zeros((n*n, 1))
    it = 0
    for var_a in v_alt[0:(n*n)]:
        #print var.VarName, var.X
        P_alt[it, 0] = var_a.X
        it += 1
    P_alt = P_alt.reshape((n, n)).transpose()
    session.putvalue('mod_Phat', P_alt)
    perm_equal = np.allclose(P, P_alt)
    print "Two permutation solutions are same? "+str(perm_equal)
    print "Nonzero elements are "+str(np.nonzero(P-P_alt))
    run_script = '[mod_tm_m'+str(mc) +', mod_delta_m'+str(mc) +'] = PermMat2TrueMatch(mod_Phat, true_perm, m)'
    session.run(run_script)
    truematch_mod.append(session.getvalue("mod_delta_m"+str(mc)))
    runtime_mod.append(gmodel_relaxed_alt.Runtime)
    #tuning = gmodel_relaxed.tune()
    #gmodel_relaxed.getTuneresult(0)
    #gmodel_relaxed.optimize()
    # optimal_runtime[mc] =  gmodel_relaxed.Runtime

tm_df = pd_df({
    "tm": truematch_orig + truematch_mod,
    "runtime": runtime_orig + runtime_mod,
    "orig_or_mod": ["orig"]*nmc + ["mod"]*nmc
})

p_val_truematch_diff = st.ttest_ind(truematch_orig, truematch_mod)
p_val_timediff = st.ttest_ind(runtime_orig, runtime_mod)


#print gg.ggplot(tm_df, aes('orig_or_mod', 'tm')) + \
#  gg.geom_line(colour='steelblue')


#print gg.ggplot(tm_df, aes('orig_or_mod', 'runtime')) + \
#  gg.geom_line(colour='steelblue')







