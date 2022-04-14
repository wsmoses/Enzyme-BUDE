#!/usr/bin/python3.8
import os

def printfun(numthreads, blocklist,itercount):
  for s in blocklist:
      for mode in ["-forward","-gradient"]:
        os.system("OMP_NUM_THREADS={}  taskset -c 0-{} numactl -i all ~/enzyme-sc22/BUDE/openmp/ompOpt-single{}.exe -n {} -i {} > omp-single{}_{}_{}_{}.txt".format(numthreads, numthreads-1, mode, s*s*s, itercount, mode,numthreads, itercount, s))

itercount=8
#Strong scaling
for numthreads in [1,8,16,24,32,40,48,56,64]:
  printfun(numthreads, [96], itercount)

