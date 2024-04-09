## Command operation after nscf calculation is complete

  It is important to note that older versions of QE 5.0 or 5.1 should be used,

```
module load espresso/5.0.2
```

  as well as the corresponding 6 cores,

```
mpirun -np 6 projwfc.x < projwfc.inp > dos.out
```

  sum_tmp.ipynb is to read the DOS data from the tmp_ID folder and summarize it.
  The summarized data is read by the code of ipynb in the Plotting folder to draw the TDOS plot.
