# gencrit
For a given choice of integer parameters 4 <= k <= n, the purpose of this program is to list all simple graphs on n vertices with no isolated vertices that are edge-k-critical, meaning that the graph has chromatic number k and the removal of any edge from the graph results in a graph with chromatic number k-1. All the outputs for n <= 13 as well as the output for n=14, k=4 can be found over at http://users.cecs.anu.edu.au/~bdm/data/graphs.html. It should be noted that unless k is close to n, a great deal of computing power would be required to compute any list that is not already available in the preceding link. For example, the computation for n=14 and k=5 would probably take something like 50 years of CPU time.

This program relies on the computational graph theory library nauty (https://users.cecs.anu.edu.au/~bdm/nauty/). The required dependencies are geng.c and nautyW1.a. To build the program, given the desired value of k, run:
```
$(CC) -o gencrit -Ofast -DMAXN=WORDSIZE -mpopcnt -march=native -DTARGET_CHI=k -DWORDSIZE=32 -DPRUNE=prune_crit -DPREPRUNE=preprune_crit geng.c gencrit.c nautyW1.a
```
To run the program on a single thread, execute the following command:
```
./gencrit -d[k-1] [n] [outputfile]
```
Or, alternatively, if you wish to use parallelism to speedup the computation, for a given choice of a natural number m and for every 0 <= i < m, you can run:
```
./gencrit -d[k-1] [n] [i]/[m] [outputfile]
```

Merging all the files obtained in this way will produce the same list as if the previous command had been run. The graphs are stored in graph6 format: refer to formats.txt in the nauty package for more information.

Many thanks to Brendan McKay for helpful guidance and for hosting the resulting files.
