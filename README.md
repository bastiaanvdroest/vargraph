# vargraph
Making variation graphs and finding virus strains in the graphs.

- graph_abundance.py: Map reads to nodes of a variation graph and calculate abundance of the nodes.
- graph_analysis.py:  Analysis of variation graph, with optimization step for finding virus strains
- graph_poa.py:       Perform partial order alignment and compress the initial variation graph
- gurobi_optimize.py: Optimization code, with the complete optimization problem
- gurobi_optmize2.py: Optimization code, with the problem split in first optimizing the binary matrix and then optimizing for f
- permutation.py:     Permute order of contigs in POA and calculate mean length of graph
