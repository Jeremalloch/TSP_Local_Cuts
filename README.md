# Using DNN to Find Local Cut Vertex Clusters
This repository contains a version of the Concorde TSP solver that has been
modified to use a deep neural network (DNN) to cluster vertices for the local
cuts routine within Concorde.  This project is being conducted as part of 
CO 759 - Deep Learning in Computational Discrete Optimization, taught by Bill
Cook at the University of Waterloo, Winter 2018.

## Concorde TSP Solver
Concorde is a computer code for the symmetric traveling salesman problem (TSP) and some related network optimization problems. The code is written in the ANSI C programming language and it is available for academic research use; for other uses, contact William Cook for licensing options.

Concorde's TSP solver has been used to obtain the optimal solutions to the full set of 110 TSPLIB instances, the largest having 85,900 cities.

The Concorde callable library includes over 700 functions permitting users to create specialized codes for TSP-like problems. All Concorde functions are thread-safe for programming in shared-memory parallel environments; the main TSP solver includes code for running over networks of UNIX workstations.
