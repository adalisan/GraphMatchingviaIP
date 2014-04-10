__author__ = 'Sancar'

from gurobipy import *
import os
import copy
import pymatlab
import numpy as np

pr_dir = os.environ.get('PROJECT_DIR')
session = pymatlab.session_factory()

#pr_dir = 'G:\Sancar\Documents\projects'
#pr_dir = '~/projects/'
print pr_dir

root_dir  = os.path.join(pr_dir ,'SeededGraphMatch')
models_dir = os.path.join(root_dir ,'optim_models')


for file in os.listdir(models_dir):
    if file.endswith(".lp"):
        gmodel = read(fname)
        gmodel.update()


        allConstr = gmodel.getConstrs( )
        relaxedConstr_22= allConstr[0:(n*n)]
        relaxedConstr_12= allConstr[(n*n):(n*n+2*m*n)]
        pen_Constr_22 = [1.0]*len(relaxedConstr_22)
        pen_Constr_12 = [1.0]*len(relaxedConstr_12)
        gmodel_relaxed = gmodel.copy()
        gmodel_relaxed.feasRelax(0, False, None, None, None, relaxedConstr_22+relaxedConstr_12, pen_Constr_22+pen_Constr_12)
        gmodel_relaxed.update()
        gmodel_relaxed.write(fname+'relaxed.lp')
        optresults = gmodel_relaxed.optimize()

        v = gmodel.getVars()
        print v
        for var in v:
            print var