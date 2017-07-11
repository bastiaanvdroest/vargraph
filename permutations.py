#!/usr/bin/python

from os import remove
from time import clock
from random import sample
from subprocess import call
from itertools import permutations
from statistics import mean, stdev
import argparse
import graph_poa
import datetime

def random_permutation(iterable, r=None):
    #"Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(sample(pool, r))

parser = argparse.ArgumentParser(prog='python permutations.py', description='Permute the order of contigs and compute the variation graph')
parser.add_argument('contigfile',type=str,help='The file containing the contigs, format is .fasta')
parser.add_argument('--r',type=str,default='results_permutations.txt',help='The file for storing the results of the permutation test, default is results_permutations.txt')
parser.add_argument('--s',type=str,default='tmp.txt',help='The file for storing the POA results, default is tmp.txt')
parser.add_argument('--m',type=str,default='blosum80.mat',help='The alignment scoring matrix, default is blosum80.mat')
parser.add_argument('--g',type=str,default='graph.gt',help='The file for storing the graph, default is graph.gt')
parser.add_argument('--p',type=str,default='graphviz.png', help='PNG file to store the visualisation of the graph, default is graphviz.png')
parser.add_argument('--make_graph',action='store_true', help='Option to make a visualization of the graph')

args = parser.parse_args()
contigfile = args.contigfile

strainline = {}
with open(contigfile) as file:
    for line in file:
        line = line.strip()
        if line.startswith('>'):
            strain = line.split('>')[1]
        else:
            strainline[strain] = line

perms = permutations(strainline.keys())

#random_picks = random_permutation(perms, 3)
length = []
n = 1
pngfile = args.p.split('.')[0]+'_'
for p in perms:
    with open('segfile'+str(n)+'.fasta','w') as file:
        for strain in p:
            file.write('>'+strain+'\n')
            file.write(strainline[strain]+'\n')
    args.contigfile = file.name
    args.p = pngfile+str(n)+'.png'
    length.append(graph_poa.main(args))
    remove(file.name)
    n += 1

with open(args.r,'w') as file:
    file.write('Results from permutations.py at %s\n\n' %now.strftime("%Y-%m-%d %H:%M"))
    file.write('Contigfile             = %s\n' %contigfile)
    file.write('Scoring matrix         = %s\n' %args.m)
    file.write('Number of permutations = %d\n' %(n-1))
    file.write('Runtime                = '+str(clock())+'\n\n')
    file.write('The results of the permutations are:\n')
    file.write('\tmean length of the graph     = %d\n' %mean(length))
    file.write('\tstd error of length of graph = %d\n' %stdev(length))
