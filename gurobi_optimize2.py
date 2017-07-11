#!/usr/bin/python

from numpy import array, c_
from gurobipy import *
from graph_tool.all import *
from tqdm import tqdm
from multiprocessing import Pool

###############################################################################
##
## constr_v(m, g, a, nstrains):
##
## Define constraints for the A matrix with binaries. In this matrix a 0 states
## that a node is not contained in a strain and a 1 that it is contained.
## A[i,j] for node i and strain j depends on the in-neighbors of i.
## If node i has no in-neighbors, then a new strain should start in this node.
## If node i has 2 or more in-neighbors, then A[i,j] is the sum of the
##  in-neighbors of i.
## If node i has 2 or more out_neighbors, then the sum of A[k,j] with k in K
##  the out_neighbors of node i, must be between A[i,j] and 1.
##
## Input:
##   - m        := the optimization model
##   - g        := the variation graph
##   - A        := the matrix with binaries for nodes and strains
##   - nstrains := the number of strains
##
###############################################################################

def constr_v(m, g, A, nstrains):
    start_vs = []
    for v in g.vertices():
        # print(v)
        # print(A[int(v),:])
        if len(list(v.in_neighbours())) == 0:  #Save nodes with no in-neighbors
            start_vs.append(int(v))
        l_in = len(list(v.in_neighbours()))   #Determine number of in-neighbors
        l_out = len(list(v.out_neighbours())) #Determine number of out_neighbors
        # print(v, l_in)

        if l_in > 1:                   #2 or more in-neighbors
            for i in list(v.in_neighbours()):
                if len(list(i.out_neighbours())) == 1:
                    for j in range(nstrains):
                        m.addConstr(A[int(v),j] >= A[int(i),j])

        elif l_in == 1:
            w = list(v.in_neighbours())[0]
            w_out = len(list(w.out_neighbours()))
            if w_out == 1:
                strains_v = g.vp.strain[v]
                strains_w = g.vp.strain[w]
                if len(strains_v) > len(strains_w):
                    for j in range(nstrains):
                        m.addConstr(A[int(v),j] >= A[int(w),j])
                elif len(strains_v) < len(strains_w):
                    for j in range(nstrains):
                        m.addConstr(A[int(v),j] <= A[int(w),j])
                else:
                    for j in range(nstrains):
                        m.addConstr(A[int(v),j] == A[int(w),j])
            else:
                npath = 0
                geq = False
                no_path = False
                for i in list(w.out_neighbours()):
                    if i.in_degree() > 1:
                        geq = True
                        for j in list(i.in_neighbours()):
                            v_list, e_list = shortest_path(g,w,j)
                            if len(e_list) == 0:
                                no_path = True
                    v_list, e_list = shortest_path(g,v,i)
                    if len(e_list) > 0:
                        npath += 1
                        gap_node = i

                if npath > 0:
                    w_out = npath+1
                else:
                    w_out = 1
                for j in range(nstrains):
                    sum_strain = 0
                    for i in list(w.out_neighbours()):
                        sum_strain += A[int(i), j]

                    # m.addConstr(sum_strain >= A[int(w),j])
                    if geq:
                        # print(sum_strain ,'>=', A[int(w),j])

                        if w_out > 1 and no_path == False:
                                m.addConstr(sum_strain <= w_out*A[int(w),j])
                        elif w_out == 1 and no_path:
                            m.addConstr(sum_strain >= w_out*A[int(w),j])
                            m.addConstr(sum_strain <= 1)
                        else:
                            m.addConstr(sum_strain >= w_out*A[int(w),j])
                    else:
                        m.addConstr(sum_strain == w_out*A[int(w),j])

    # for i in range(len(start_vs)):    #Nodes with no in-neighbors are start of
    #     for j in range(start_vs[i]):  #a new strain
    #         m.addConstr(A[j,i] == 0)
    #     m.addConstr(A[start_vs[i],i] == 1)

###############################################################################
##
## f, a, objVal = optimize(b, nvert, nstrains, g, options):
##
## Optimization function calculates the frequencies of the strains and
## determines which node is contained in which strain.
## It also calulcates the number of strains that should be used to minimize
## the difference between the abundance of a node and the sum of frequencies
## of all strains in that node.
##
## The function needs as input:
##   b        := the abundances of the nodes
##   nvert    := the number of nodes
##   nstrains := the number of strains
##   g        := the variation graph
##
## The function has as output:
##   f        := the frequencies of the strains
##   a        := the matrix determining whether a node is included in a strains
##   objVal   := the minimal value of the objective function
##
###############################################################################

def optimize(model_file, b, nvert, nstrains, g, options):
    M = 1000000

    # Read in model
    m = read(model_file)
    mvars = list(m.getVars())
    e1 = []
    e2 = []
    for i in range(len(mvars)):
        if i < nvert:
            e1.append(mvars[i])
        else:
            e2.append(mvars[i])

    # Create variables
    alpha = m.addVars(list(range(nvert*nstrains)), vtype=GRB.BINARY, name = 'alpha')

    m.update()

    A = []                                  #Put the alpha variables in
                                            #matrix form
    for i in range(nvert*nstrains):
        A.append(alpha[i])
    A = array([A]).reshape(nvert,nstrains)

    # Set objective:
    obj = LinExpr()
    for i in range(nvert):             #Take the sum over all nodes
        obj += e2[i]/max(1,b[i,0])
    m.setObjective(obj, GRB.MINIMIZE)

    m.update()

    # Add constraint:
    for i in range(nvert):
        #Sum of the frequencies of strains is the same as the sum of the
        #contigs in a node
        m.addConstr(sum(A[i,:]) + e1[i] == b[i,0])

        #Every node must be in at least one strain
        m.addConstr(sum(A[i,:]) >= 1)

    # tot_freq = 0
    # for i in range(nstrains):
    #     tot_freq += f[i]
    # m.addConstr(tot_freq >= 1)
    constr_v(m, g, A, nstrains)   # Define the constraints for the A matrix

    m.update()
    m.Params.LogToConsole = 0
    m.Params.Threads = 0
    m.Params.NumericFocus = 0
    m.Params.PoolSearchMode = 0
    m.Params.PoolSolutions = 10
    #Minimize the model for the given objective function and constraints
    m.optimize()
    # m.computeIIS()
    # m.write('model.ilp')

    if m.status == GRB.Status.OPTIMAL:
        A = []
        for v in m.getVars():
            if 'alpha' in v.varName:
                A.append(int(v.xn))
        A = array(A).reshape(nvert,nstrains)

        del m

        m = read(model_file)
        mvars = list(m.getVars())
        e1 = []
        e2 = []
        for i in range(len(mvars)):
            if i < nvert:
                e1.append(mvars[i])
            else:
                e2.append(mvars[i])

        f = m.addVars(list(range(nstrains)), lb=0, vtype=GRB.CONTINUOUS, name = 'f')
        m.update()

        obj = LinExpr()
        for i in range(nvert):             #Take the sum over all nodes
            obj += e2[i]/max(1,b[i,0])
        m.setObjective(obj, GRB.MINIMIZE)
        m.update()
        # Add constraint:
        for i in range(nvert):
            sum_vert = 0
        #Sum of the frequencies of strains is the same as the sum of the
        #contigs in a node
            for j in range(nstrains):
                sum_vert += A[i,j]*f[j]
            m.addConstr(sum_vert + e1[i] == b[i,0])

        m.update()
        m.Params.LogToConsole = 0
        m.Params.Threads = 0
        m.Params.NumericFocus = 0
        m.Params.PoolSearchMode = 0
        m.Params.PoolSolutions = 10
        #Minimize the model for the given objective function and constraints
        m.optimize()

        if m.status == GRB.Status.OPTIMAL:
            for i in reversed(range(m.solcount)):
                #Create lists with the frequencies optimized by the model
                f = []

                if m.solcount > 1:
                    m.setParam('SolutionNumber',i)
                    objVal = m.PoolobjVal

                objVal = m.objVal
                for v in m.getVars():
                    if 'f' in v.varName:
                        f.append(v.X)

                if options['display']:      #If the display option is on, print
                    print('\nStatus: optimal')      # - status = optimal
                    print('\nSolution number %d' %(i+1))
                    print('\nObjective value: %g' % objVal) # - the obj value

                    print('\nFrequencies of the strains are:') # - the frequencies
                    for v in m.getVars():
                        if 'f' in v.varName:
                            print(' %s = %g' % (v.varName, v.x))

                    print('\nV matrix is:') # - the v matrix with the binaries
                    print(A,'\n')

                    seq = {}
                    for j in range(len(A[0,:])):
                        seq[j] = ''
                        for k in range(len(A[:,0])):
                            if int(A[k,j]) == 1:
                                seq[j] += g.vp.seq[k]
                    print('The sequences of the strains are:\n')
                    for k in range(len(A[0,:])):
                        print(' s[%d] ='%k,seq[k],'\n')

            return(f, A, objVal)
        else:
            return([], [], 1e9) #If the solution is not optimal, return empty lists
                                #for f and a and set obj value at 100.
    else:
        return([], [], 1e9)
    del m

if __name__ == '__main__':
    optimize(b, nvert, nstrains, g, options)
