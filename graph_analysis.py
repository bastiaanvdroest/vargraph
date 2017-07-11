#!/usr/bin/python

import sys
import os.path
from graph_tool.all import *
from numpy import *
from gurobipy import *
from gurobi_optimize2 import optimize
import time
import multiprocessing
from tqdm import tqdm
import argparse
import datetime

def give_abun(args, g, nvert, ncontigs):
    if hasattr(args,'a'):
        vprop = g.new_vertex_property('double')
        g.vp.abun = vprop
        with open(args.a) as file:
            data = file.readlines()
        for i in range(len(data)):
            data[i] = int(data[i])
        b = []

        b_neg = []
        vfilt = g.new_vertex_property('bool')
        for i in range(nvert):
            if float(data[i]) < 20000*args.f:
                vfilt[i] = False
                b_neg.append(i)
            else:
                vfilt[i] = True
            b.append(float(data[i]))
            g.vp.abun[i] = data[i]
        del data
        g = GraphView(g, vfilt, skip_vfilt=True)
        g.purge_vertices()
        nvert = len(list(g.vertices()))
        b = array([delete(b,b_neg).tolist()]).reshape(nvert,1)

    else:
        abun_contig = input('Give the list of abundances of the contigs:').split(',')
        args.a = abun_contig
        b = zeros((nvert,1))
        for v in g.vertices():
            for j in range(ncontigs):
                if str(j) in g.vp.strain[v]:
                    b[int(v),0] += float(abun_contig[j])
    return b, g, nvert, args

def give_results(objVal,f,V,g,args):
    now = datetime.datetime.now()
    with open(args.r,'w') as file:
        file.write('Results from graph_analysis.py at %s\n\n' %now.strftime("%Y-%m-%d %H:%M"))
        file.write('Input is:\n\n')
        file.write('  Variation graph = %s\n' %args.graphfile)
        file.write('  Abundances      = %s\n\n' %args.a)
        file.write('Output is:\n\n')
        file.write('  Runtime         = %g\n' %time.clock())
        file.write('  Status          = optimal\n')
        file.write('  Objective value = %g\n' %objVal)
        print('\nObjective value: %g' % objVal) # - the obj value

        file.write('  Frequencies of the strains are:\n')
        print('\nFrequencies of the strains are:') # - the frequencies
        for i in range(len(f)):
            print('f[%d] = %g'%(i,f[i]))
            file.write('    f[%d] = %g\n'%(i,f[i]))

        file.write('\n  V matrix is:\n')
        print('\nV matrix is:') # - the v matrix with the binaries
        for slice in V:
            slice = str(slice.tolist())
            file.write('    %s\n'%slice)
        print(V,'\n')

        seq = {}
        for j in range(len(V[0,:])):
            seq[j] = ''
            for i in range(len(V[:,0])):
                if int(V[i,j]) == 1:
                    seq[j] += g.vp.seq[i]
        print('The sequences of the strain are:\n')
        file.write('\n  Sequences of the strains are:\n')
        for i in range(len(V[0,:])):
            if f[i] > 0:
                print(' s[%d] ='%i,seq[i],'\n')
                file.write('    s[%d] = %s\n' %(i,seq[i]))
            for j in range(len(V[0,:])):
                same = True
                for k in range(len(V[:,0])):
                    if V[k,i] != V[k,j]:
                        same = False
                if same:
                    print('strains are the same: ',(i,j))

def opt_step(b,nvert,ncontigs,g,options):
    #Calculate the supporting points
    k = []
    P = []
    print('\nCalculating the supporting points for approximating the problem:')
    for i in tqdm(range(len(b))):
        P.append([])
        k.append(ceil(4*b[i,0]/sqrt(0.001*max(1,b[i,0]))))
        for j in range(int(k[i])):
            P[i].append(-2*b[i,0]+(2*j-1)*sqrt(0.001*max(1,b[i,0])))
    del k

    # Run the optimalizations
    # pool = multiprocessing.Pool(processes = 1)
    #
    # tasks = []
    # results = []
    m = Model('lp')
    e1 = m.addVars(list(range(nvert)), vtype=GRB.CONTINUOUS, name='e1')
    e2 = m.addVars(list(range(nvert)), lb=0, vtype=GRB.CONTINUOUS, name='e2')
    m.update()

    print('\nSetting constraints for the supporting points:')
    for i in tqdm(range(nvert)):
        for j in range(len(P[i])):
            # print(2*P[i][j]*e1[i],'-',e2[i],'<=',P[i][j]**2)
            m.addConstr(2*P[i][j]*e1[i] - e2[i] <= P[i][j]**2)
    m.update()
    model_file = 'model.mps'
    m.write(model_file)
    del m

    min_objVal = 1e9

    print('\nPerforming the optimalization steps: ')
    for nstrains in tqdm(range(1,ncontigs+1)):
        f, V, objVal = optimize(model_file, b, nvert, nstrains, g, options)
    # for nstrains in range(1,ncontigs+1):
    #     tasks.append((b, k, p, nvert, nstrains, g, options))
    #
    # for t in tasks:
    #     results.append(pool.apply_async(optimize, t))
    #
    # for result in tqdm(results):
    #     f, V, objVal = result.get()
    #     objvs.append(objVal)
        if objVal < min_objVal:
            min_objVal = objVal
            min_f = f
            min_V = V
        del f
        del V
    return min_objVal, min_f, min_V

def visual_graph(g):
    gprops = {'K': 0.5, 'ordering': 'out', 'rank': 'source', 'rankdir': 'LR', 'ratio': 'auto'}
    vprops_seq = {'label': g.vertex_index, 'shape': 'box', 'fillcolor': g.vp.col, }#'xlabel': g.vp.abun}
    eprops = {'arrowsize': 2.0, 'dir': 'forward'}

    graphviz_draw(g, output=args.p, layout = 'dot', size = (500,500), vsize = 1, overlap=False, gprops = gprops, vprops= vprops_seq, eprops = eprops)

def main(args):
    options = {'optimize':args.opt, 'make_graph':args.make_graph, 'display':args.d}
    if hasattr(args,'r'):
        options['results'] = True
    else:
        options['results'] = False

    g = load_graph(args.graphfile)
    vfilt = g.new_vertex_property('bool')
    # for i in range(0,50):
    #     vfilt[i] = True
    # g = GraphView(g,vfilt)
    nvert = len(list(g.vertices()))
    ncontigs = 0
    with open(args.contigfile) as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                ncontigs += 1

    b, g, nvert, args = give_abun(args, g, nvert, ncontigs)

    if options['make_graph']:
        visual_graph(g)

    if options['optimize']:
        objVal, f, V = opt_step(b,nvert,ncontigs,g,options)
        if options['results']:
            give_results(objVal,f,V,g,args)

    print('Runtime is', time.clock())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='python graph_analysis.py', description='Find abundances and sequences of strains in a given variation graph')
    parser.add_argument('contigfile',type=str,help='The file containing the contigs, format is .fasta')
    parser.add_argument('graphfile',type=str,help='The file containing the variation graph, format is .gt')
    parser.add_argument('-a',type=str,default=argparse.SUPPRESS,help='The file with the abundances of the nodes, format is .txt')
    parser.add_argument('-opt',action='store_true',help='Perform the optimization to find strains')
    parser.add_argument('-r',nargs='?',const='results_optimization.txt',default=argparse.SUPPRESS,help='Show and store the results of the optimization')
    parser.add_argument('-d',action='store_true',help='Show the intermediate results')
    parser.add_argument('-p',type=str,default='graphviz.png', help='PNG file to store the visualisation of the graph, default is graphviz.png')
    parser.add_argument('-f',nargs='?',const=0.1,default=0,help='Filter out the nodes of the variation graph that are below the given percentage F of the read coverage')
    parser.add_argument('-make_graph',action='store_true', help='Option to make a visualization of the graph')

    args = parser.parse_args()
    main(args)
