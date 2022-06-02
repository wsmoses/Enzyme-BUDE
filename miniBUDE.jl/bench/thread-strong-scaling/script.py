#!/usr/bin/python3.8
import os

def printfun(numthreads, itercount):
      os.chdir("../../")
      for mode in ["false","true"]:
        os.system("julia --project=Threaded --threads={} src/Threaded.jl --enzyme {} --iterations {}  > omp-single-{}_{}_{}.txt".format(numthreads,mode,itercount,mode,numthreads,itercount))
        os.system(" mv *.txt bench/thread-strong-scaling/")
      os.chdir("bench/thread-strong-scaling")

itercount=8
#Strong scaling
for numthreads in [1,8,16,24,32,40,48,56,64]:
  printfun(numthreads, itercount)
