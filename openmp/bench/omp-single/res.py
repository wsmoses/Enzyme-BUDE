#!/usr/bin/python3.8
import os

def printfun(numthreads, blocklist,itercount,mode):
     for s in blocklist:
        os.system("echo -n {},{}, >> results.txt".format(numthreads,s))
        os.system("grep  Average omp-single{}_{}_{}_{}.txt >> results.txt".format(mode,numthreads, itercount, s))
     os.system("sed -i \"s/- Average time:   //g\" results.txt")
     os.system("sed -i \"s/ ms//g\" results.txt")

itercount=100
os.system("rm -f results.txt")
#Strong scaling
for mode in ["-forward","-gradient"]:
  for numthreads in [1,8,16,24,32,40,48,56,64]:
    printfun(numthreads, [96], itercount,mode)

