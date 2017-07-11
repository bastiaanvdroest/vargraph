#!/usr/bin/python

import sys
import os.path
from graph_tool.all import *
import time
import math
from subprocess import call
import argparse


#####################################################################
## Function make_graph(contigfile, storefile, scoring_matrix):
##
## This function makes a graph with only one nucleotide as the
## sequence of each node. To make this graph, the function performs
## a partial order alignment on the given contigs and the blosum80
## scoring matrix. It then translates the results into a graph using
## the graph_tool library.
##
## Input for the function is:
##    - contigfile:        File with all the contigs in .fasta format
##    - storefile:         Empty file where the results can be stored
##    - scoring_matrix:    The alignment scoring matrix, default is
##                         the blosum80 matrix
##
## Output of the function is:
##    - g:                 The variation graph made with the given
##                         contigs, composed out of nodes with one
##                         nucleotide
##    - v_start:           List of vertices where the contigs start
#####################################################################
def make_graph(contigfile, storefile, scoring_matrix = 'blosum80.mat'):
    call(' '.join(['poa -read_fasta', contigfile, '-po', storefile, scoring_matrix]), shell=True)

    g = Graph(directed=True)
    vprop = g.new_vertex_property('string')
    g.vp.seq = vprop
    vprop = g.new_vertex_property('vector<string>')
    g.vp.strain = vprop

    f = open(storefile)

    c_start = []
    v_start = []
    c=0
    for line in f:
        line = line.strip()
        base = line.split(':')[0]
        if base in ['A','T','G','C','N']:
            c+=1
            rest = line.split(':')[1]
            v = g.add_vertex()
            g.vp.seq[v] = base
            for i in range(len(rest)):
                if rest[i] == 'L':
                    node = ''
                    n = i+1
                    while rest[n] not in ['L','S','A'] and n < len(rest)-1:
                        node += rest[n]
                        n += 1
                    e = g.add_edge(node,v)
                elif rest[i] == 'S':
                    seq = ''
                    n = i+1
                    while n <= len(rest)-1:
                        if rest[n] not in ['L','S','A']:
                            seq += rest[n]
                            n += 1
                        else:
                            break
                    if seq not in c_start:
                        if int(v) not in v_start:
                            v_start.append(int(v))
                        c_start.append(seq)
                    g.vp.strain[v].append(seq)

    return (g,v_start)

def compress_graph(mg_output):
    g, v_start = mg_output
    vprop = g.new_vertex_property('string')
    g.vp.col = vprop

    n = 0
    while n < len(list(g.vertices()))-1:
        v = list(g.vertices())[n]
        g.vp.col[v] = 'red'

        if v.out_degree() == 0:
            if v.in_degree() == 0:
                for i in range(len(v_start)):
                    if v_start[i] >= int(v):
                        v_start[i] += -1
                g.remove_vertex(v)
            else:
                n += 1
        elif v.out_degree() == 1:
            w = list(v.out_neighbours())[0]
            if w.in_degree() == 1:
                if int(w) not in v_start:
                    strains_v = list(g.vp.strain[v])
                    strains_w = list(g.vp.strain[w])
                    if len(strains_v) == len(strains_w):
                        g.vp.seq[w] = ''.join([g.vp.seq[v],g.vp.seq[w]])
                        g.vp.strain[w] = list(strains_v)
                        g.vp.col[w] = g.vp.col[v]
                        v_in_neighbours = list(v.in_neighbours())
                        while v_in_neighbours != []:
                            e = g.add_edge(v_in_neighbours[0],w)
                            del(v_in_neighbours[0])
                        for i in range(len(v_start)):
                            if v_start[i] > int(v):
                                v_start[i] += -1
                        g.remove_vertex(v)
                    else:
                        n += 1
                else:
                    n += 1
            else:
                n += 1
        else:
            n += 1

    length_g = 0
    for v in g.vertices():
        length_g += len(g.vp.seq[v])
        if v.out_degree() == 0 or v.in_degree() == 0:
            g.vp.col[v] = 'cyan'
        if v in v_start:
            g.vp.col[v] = 'cyan'

    print('The total sequence length of the graph is',length_g)
    return (g, length_g)

def main(args):
    make_graph_output = make_graph(args.contigfile, args.s, args.m)
    g, length_g = compress_graph(make_graph_output)
    g.save(args.g)

    gprops = {'K': 0.5, 'ordering': 'out', 'rank': 'source', 'rankdir': 'LR', 'ratio': 'auto'}
    vprops_seq = {'label': g.vp.seq, 'xlabel':g.vp.strain,'shape': 'box', 'fillcolor': g.vp.col}
    # vprops_strain = {'label': g.vp.strain, 'shape': 'box', 'fillcolor': g.vp.col}
    eprops = {'arrowsize': 2.0, 'dir': 'forward'}

    # vfilt = g.new_vertex_property('bool')
    # for i in range(0,14):
    #     vfilt[i] = True
    # u = GraphView(g,vfilt)

    if args.make_graph:
        graphviz_draw(u, output= args.p, layout = 'dot', size = (500,500), vsize = 1, overlap=False, gprops = gprops, vprops= vprops_seq, eprops = eprops)
    #
    # graphviz_draw(g, output='graphviz_strain.png', layout = 'dot', size = (500,500), vsize = 1, overlap=False, gprops = gprops, vprops= vprops_strain, eprops = eprops)

    print('Runtime is', time.clock())

    return length_g

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='python graph_poa.py', description='Perform partial order alignment and compress the variation graph')
    parser.add_argument('contigfile',type=str,help='The file containing the contigs, format is .fasta')
    parser.add_argument('--s',type=str,default='tmp.txt',help='The file for storing the POA results, default is tmp.txt')
    parser.add_argument('--m',type=str,default='blosum80.mat',help='The alignment scoring matrix, default is blosum80.mat')
    parser.add_argument('--g',type=str,default='graph.gt',help='The file for storing the graph, default is graph.gt')
    parser.add_argument('--p',type=str,default='graphviz.png', help='PNG file to store the visualisation of the graph, default is graphviz.png')
    parser.add_argument('--make_graph',action='store_true',help='Option to visualize the graph')

    args = parser.parse_args()

    main(args)
