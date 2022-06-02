#!/usr/bin/python3.8
import os
import pathlib
scriptdir = pathlib.Path(__file__).parent.resolve()

def printfun(numthreads, blocklist,itercount):
  for s in blocklist:
      for mode in ["-forward","-gradient"]:
        os.system("OMP_NUM_THREADS={}  taskset -c 0-{} numactl -i all ".format(numthreads, numthreads-1) + str(scriptdir) + "/../../ompOpt-single{}.exe -n {} -i {} --deck {} > ".format(mode, s*s*s, itercount, str(scriptdir)+"/../data/bm1/")+str(scriptdir)+"/omp-single{}_{}_{}_{}.txt".format(mode,numthreads, itercount, s))

itercount=100
#Strong scaling
for numthreads in [1,8,16,24,32,40,48,56,64]:
  printfun(numthreads, [96], itercount)
