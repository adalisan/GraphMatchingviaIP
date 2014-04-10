import igraph

from igraph import Graph as gr
from gurobipy import *
import numpy  as np
import copy as cp

#graph_1 = gr.Read_GML("graph.pair.1.gml")
#graph_2 = gr.Read_GML("graph.pair.2.gml")

graph_1 = gr.Erdos_Renyi(n=40, p=0.5)
A_Matrix = graph_1.get_adjacency()
A_dat= A_Matrix.data
A = np.array(A_dat)
#B = graph_1.get_adjacency()
B_matrix = cp.deepcopy(A)
#B= np.array(B_matrix)
B= B_matrix
m = 10


mplusn = graph_1.vcount()
n = mplusn-m
seed_indices = range(m)
nonseed_indices = range(m,mplusn)

trueperm = cp.deepcopy(nonseed_indices)
#np.random.shuffle(trueperm)
shuff_indices = seed_indices+trueperm
B= B[shuff_indices, :][:, shuff_indices]

A11=A[seed_indices, :][:, seed_indices]
A12=A[seed_indices, :][:, nonseed_indices]
A21=A[nonseed_indices, :][:, seed_indices]
A22=A[nonseed_indices, :][:, nonseed_indices]
B11=B[seed_indices, :][:, seed_indices]
B12=B[seed_indices, :][:, nonseed_indices]
B21=B[nonseed_indices, :][:, seed_indices]
B22=B[nonseed_indices, :][:, nonseed_indices]
#
# # does numpy store in column order or row order
vecB12=B12.reshape(m*n,1)
vecA21=A21.reshape(m*n,1)
#
# ident_n = np.identity(n,float)
# ident_nsq = np.identity(n^2+2*m*n, float)
# M11=vstack( kron(ident_n,A22)-kron(B22.transpose(),ident_n),
#       kron(ident_n,A12), kron(B21.transpose(),ident_n) )
# M21 = vstack( kron(ident_n,ones((1,n),float) , kron(ones((1,n),float),ident_n)))
# M = vstack( hstack( M11, ident_nsq,  -ident_nsq ),
#             hstack( M21 , zeros((2*n,2*n^2+4*m*n),float)  ) )
# f = hstack( zeros(1,n^2), ones((1,2*n^2+4*m*n) ),float)
# Mineq= zeros((0,3*n^2+4*m*n),float)
# b=vstack(zeros(n^2,1) , vecB12 , vecA21 , ones((2*n,1)),float)
#
# sense_eq = tile('=',n^2+2*m*n+2*n, 1)
# sense_ineq = repmat('<',0, 1);
# sense = [sense_eq; sense_ineq];
# %Binary and Real Variables
# vtype1 = repmat('B',n^2, 1);
# vtype2 = repmat('C',2*n^2+4*m*n,1);
# vtype = [vtype1; vtype2];
#
# sgm_solver_mod = Model("SGM")
# sgm_solver_mod.se A = M;
# model.obj = f;
# model.modelsense = 'min';
# model.rhs = b;
# model.sense = sense;
# model.vtype = vtype;




# Create a new model
gmodel = Model("SGM")
perm_mat_entries = {}
# Create variables
for i in range(n):
    for j in range(n):
        perm_mat_entries[i,j] = gmodel.addVar(vtype=GRB.BINARY,lb=0, ub=1, obj=0.0, name='P_{'+str(i)+','+str(j)+'}' )

slacks_pos = []
slacks_neg = []
for i in range(n**2+(2*m*n)):
    slacks_pos.append( gmodel.addVar(vtype=GRB.CONTINUOUS, lb=0.0, obj=1.0, name='slack_pos'+str(i)))

for i in range(n**2+(2*m*n)):
    slacks_neg.append( gmodel.addVar(vtype=GRB.CONTINUOUS, lb=0.0, obj=1.0, name='slack_neg'+str(i)))


gmodel.setAttr("ModelSense", GRB.MINIMIZE)
# Integrate new variables
gmodel.update()
# Set objective

# Add perm mat constraints: defined by lb and ub already



# Add perm_mat*slack vars constraints

# For (a*n+b)th constraint (which involves the (a*n+b)th row of M_11 = I \kron A22 -B22 \kron I
for a in range(n):
    for b in range(n):
        gmodel.addConstr(
            quicksum(perm_mat_entries[d, a]*a_22  for d in range(n) for a_22 in A22[range(n), d]) -
            quicksum(b_22*perm_mat_entries[b, c]  for c in range(n) for b_22 in B22[range(n), c]) +
            slacks_pos[a*n+b] - slacks_neg[a*n+b] == 0.0,
            name='perm_mat_diag_blk_%s_%s' % (a, b))

#First mn slacks   I \kron A_12
for a in range(n):
    for b in range(m):
        gmodel.addConstr((quicksum(perm_mat_entries[d, a]*A12[b, d]  for d in range(n))  +
                slacks_pos[(n**2)+a*m+b] - slacks_neg[(n**2)+a*m+b]) == vecB12[a*m+b], name='a_12_to_b_12_%s_%s' % (a, b))



#Second mn slacks  I \kron B_21^T
for a in range(m):
        for b in range(n):
            gmodel.addConstr(quicksum(perm_mat_entries[b, c]*B21[c, a] for c in range(n) ) +
                slacks_pos[(n**2)+(m*n)+a*n+b] - slacks_neg[(n**2)+(m*n)+a*n+b] == vecA21[a*n+b], name='b_21_to_a_21_%s_%s' % (a, b))

#questionable should slack and veca21 be a*n+b or a*m+b

for d in range(n):
        gmodel.addConstr(quicksum(perm_mat_entries[d, c] for c in range(n)) == 1.0, name='perm_mat_dbl_stoc_row_%s' % d)

for c in range(n):
        gmodel.addConstr(quicksum(perm_mat_entries[d, c] for d in range(n)) == 1.0, name='perm_mat_dbl_stoc_clmn_%s' % c)


# Add slack vars constraints
#  already defined by  lb ub of slacks

gmodel.update()

tuning= gmodel.tune()
#result = gmodel.optimize()
#
# x = result.x;
# x=round(x);
#
# P=zeros(n,n);
# for i=1:n
#     P(:,i)=x( (i-1)*n+1:i*n );
# end
# temp=P*[1:n]';
# alignment=[ [1:m] temp'+m ];

