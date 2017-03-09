# Big_GA
functions for optimizing induced sub-matrices that meet a certain fitness criterion over large matrices

Sometimes we are interested in finding a certain sub matrix within a bigger matrix that satisfies a certain 
property that we are interested in. This problem is a little hard when you consider that the sub-matrix is 
a permutation of the rows and columns of the bigger matrix and the number of ways you can choose this grows 
combinatorially. A heuristic algorithm might help in finding satisfactory solutions (may not be the best)
and I have implemented a Genetic Algorithm in a quick and dirty fashion (and possibly highly inefficient fashion)
to optimize the fitness function
