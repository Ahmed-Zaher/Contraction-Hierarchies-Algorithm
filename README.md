# Contraction-Hierarchies-Algorithm
My implementation of contraction heirarchies algorithm.

I used node level, contracted neighbours as heuristics, as well as an approximation of edge difference by assuming that we always add shortcuts.

The code was tested on a [road network graph of New York](http://www.dis.uniroma1.it/challenge9/download.shtml) (has about 250K nodes, 700K edges) and takes on average 28s for preprocessing and average query time of 3ms.
