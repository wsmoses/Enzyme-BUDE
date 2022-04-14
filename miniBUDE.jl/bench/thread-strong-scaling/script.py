#!/usr/bin/python3.8
import os

def printfun(numthreads, blocklist,itercount):
  for s in blocklist:
      for mode in ["false","true"]:
        os.chdir("/home/ubuntu/enzyme-sc22/BUDE/miniBUDE.jl/")
        os.system("julia --project=Threaded --threads={} src/Threaded.jl --enzyme {} --iterations {} --n {}  > omp-single-{}_{}_{}_{}.txt".format(numthreads,mode,itercount,s,mode,numthreads,itercount,s)) 
        os.system(" mv *.txt /home/ubuntu/enzyme-sc22/BUDE/miniBUDE.jl/bench/thread-strong-scaling/")

itercount=8
#Strong scaling
for numthreads in [1,8,16,24,32,40,48,56,64]:
  printfun(numthreads, [96], itercount)

