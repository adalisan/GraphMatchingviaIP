__author__ = 'Sancar'


from gurobipy import *
import os

pr_dir = os.environ.get('PROJECT_DIR')
pr_dir = "G:\Sancar\Documents\projects"
print pr_dir
gmodel = read(os.path.join(pr_dir ,'SeededGraphMatch','SGM_IP_prob40vt_corr_pt7_10seeds_sparse_gr.lp'))
gmodel.update()

gmodel.optimize()

#tuning= gmodel.tune()
