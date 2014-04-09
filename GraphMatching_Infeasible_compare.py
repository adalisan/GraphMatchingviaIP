#!/usr/bin/python

import sys

__author__ = 'Sancar'

from gurobipy import *
import os
import copy
import pymatlab
import numpy as np

pr_dir = os.environ.get('PROJECT_DIR')
session = pymatlab.session_factory()

#pr_dir = 'G:\Sancar\Documents\projects'
pr_dir = '~/projects/'
print pr_dir

root_dir = os.path.join(pr_dir,'SeededGraphMatch')



n = 45
m = 5

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
    gmodel_relaxed.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12,
                             pen_Constr_22+pen_Constr_12)
    gmodel_relaxed.update()
    gmodel_relaxed.write(infeas_model_fname + 'relaxed.lp')
    optresults = gmodel_relaxed.optimize()

    v = gmodel_relaxed.getVars()

    P = np.zeros((n*n, 1))
    it = 0
    for var in v[0:(n*n)]:
        #print var.VarName, var.X
        P[it,0] = var.X
        it += 1

    P = P.reshape((n,n)).transpose()

    session.putvalue('Phat', P)
    #run_script = "Psoln = DS2Perm(Phat)"
    #session.run(run_script)
    run_script = "[tm_m, delta_m] = PermMat2TrueMatch(Phat, true_perm, m)"
    session.run(run_script)
    truematch_orig.append(session.getvalue("delta_m"))
    runtime_orig.append(gmodel_relaxed.Runtime)

    pen_Constr_12 = [5.0]*len(relaxedConstr_12)

    gmodel_relaxed_alt = read(os.path.join(returndir, 'src', infeas_model_fname))
    gmodel_relaxed_alt.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, pen_Constr_22+pen_Constr_12)
    gmodel_relaxed_alt.update()


    optresults = gmodel_relaxed_alt.optimize()

    v = gmodel_relaxed_alt.getVars()

    P_alt = np.zeros((n*n, 1))
    it = 0
    for var in v[0:(n*n)]:
        #print var.VarName, var.X
        P_alt[it, 0] = var.X
        it += 1

    P_alt = P_alt.reshape((n, n)).transpose()

    session.putvalue('Phat', P_alt)
    #run_script = "Psoln = DS2Perm(Phat)"
    #session.run(run_script)
    run_script = '[tm_m, delta_m] = PermMat2TrueMatch(Phat, true_perm, m)'
    session.run(run_script)
    truematch_mod.append(session.getvalue("delta_m"))
    runtime_mod.append(gmodel_relaxed.Runtime)

    #tuning = gmodel_relaxed.tune()
    #gmodel_relaxed.getTuneresult(0)
    #gmodel_relaxed.optimize()
    # optimal_runtime[mc] =  gmodel_relaxed.Runtime


