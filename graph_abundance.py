#!/usr/bin/python

import sys
from pysam import AlignmentFile
from graph_tool.all import *
from subprocess import call
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(prog='python graph_abundance.py', description='Calculate the abundances of the nodes by mapping the reads to the contigs and count the mapped reads for each node. Abundances of the nodes will be in file abundance_[contigfile].txt')
parser.add_argument('graphfile',type=str,help='The file containing the variation graph, format is .gt')
parser.add_argument('contigfile',type=str,help='The file containing the contigs, format is .fasta')
parser.add_argument('r1',type=str,help='The file containing the reads, format is .fastq')
parser.add_argument('-r2',type=str,default=argparse.SUPPRESS,help='The file containing the other parts of the pairs for pair-end reads')
args = parser.parse_args()

g = load_graph(args.graphfile)

# Name every contig from 0 to number of contigs and make index for the contig file
name = args.contigfile.split('.')[0]

ncontigs = 0
with open(args.contigfile, 'r') as file:
    data = file.readlines()

for i in range(len(data)):
    if data[i].startswith('>'):
        data[i] = '>%s\n' % str(ncontigs)
        ncontigs += 1

with open(args.contigfile, 'w') as file:
    file.writelines(data)

nvert = len(list(g.vertices()))
print('\nNumber of contigs are: %d' % ncontigs)
print('\nNumber of vertices are: %d\n' % nvert)

call('bwa index %s' % args.contigfile, shell=True)

# Map reads against the contigs
if hasattr(args,'readfile2'):
    call('bwa mem %s %s %s | samtools view -Sb -> BAM/%s.bam' % (args.contigfile, args.r1, args.r2, name), shell=True)
else:
    call('bwa mem %s %s | samtools view -Sb -> BAM/%s.bam' % (args.contigfile, args.r1, name), shell=True)

# Sort the mapped reads and index them
call('samtools sort -T /tmp/%s.sorted -o BAM/%s.sorted.bam BAM/%s.bam' % (name, name, name),shell=True)
call('samtools index BAM/%s.sorted.bam' % name, shell=True)

bamfile = 'BAM/%s.sorted.bam' % name  # sorted.bam file with index
samfile = AlignmentFile(bamfile, 'rb')

counters = [0]*ncontigs
abundance = [0]*nvert
print('\nCounting the reads for each node:')
for v in tqdm(list(g.vertices())):
    # print(v)
    seq_length = len(g.vp.seq[v])
    for c in g.vp.strain[v]:
        counters[int(c)] += seq_length
        interval = (counters[int(c)] - seq_length, counters[int(c)])

        reads = samfile.count(c, start=interval[0], end=interval[1])
        abundance[int(v)] += reads

    abundance[int(v)] = abundance[int(v)]

file = open('abundances_%s.txt' % name,'w')
for i in range(nvert):
    file.write(str(abundance[i])+'\n')
