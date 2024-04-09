# MD calcualtion

## 1. 2x1x1, 2x2x1 and 3x3x1 super cell MD type calculation
 For the number of nodes and cores available to specify the computational 
 resources required for MD computing for different systems.
 

 Like, the cores of 24 come from 2 nodes with 12 cores per node.
 4 input files, each using the computational resources of 6 cores.

 mpirun -np 24  fdepw.x  -ni 4 -in pentacene-211

 
