## MD calcualtion

   Folders 2x1x1, 2x2x1, 3x3x1, three different system folders contain different fragment input files.

   - pentacene-211_*.in" input files describing different fragments the 2x1x1 system in the eQE calculations
   - pentacene-221_*.in" input files describing different fragments the 2x2x1 system in the eQE calculations
   - pentacene-331_*.in" input files describing different fragments the 3x3x1 system in the eQE calculations

## Step1:

 2x1x1, 2x2x1 and 3x3x1 super cell MD type calculation

 For the number of nodes and cores available to specify the computational
 resources required for MD computing for different systems.

 Like, the slurm submit file in folder 2x1x1,
 the cores of 24 come from 2 nodes with 12 cores per node.
 4 input files, each using the computational resources of 6 cores.

```
 mpirun -np 24  fdepw.x  -ni 4 -in pentacene-211
```



 
