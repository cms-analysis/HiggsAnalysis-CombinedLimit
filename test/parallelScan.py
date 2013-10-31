#!/usr/bin/env python

##  run combine -M MultiDimFit scans in multiple local jobs in parallel
##  usage: just replace "combine" with "parallelScan.py" in the command line
##         and add a "-j" to select how many threads to use
import os, sys, subprocess;
from math import ceil
from re import match

if len(sys.argv) < 3:
    print "usage: parallelScan.py <arguments to combine>  [ -j processes ]"
    exit()

jobs, points, name, method, mass, hadd = 0, 0, "Test", None, 120, False
## take out the -n and -j; find the --points
args = []
i = 1;
if sys.argv[i] == "combine": i+= 1
while i < len(sys.argv):
    if sys.argv[i] == "-j":
        jobs = int(sys.argv[i+1])
        i += 2
    elif sys.argv[i] == "-n":
        name = sys.argv[i+1]
        i += 2
    elif sys.argv[i] == "--points":
        points = int(sys.argv[i+1])
        args.append(sys.argv[i]); 
        args.append(sys.argv[i+1]); 
        i += 2
    elif "--points=" in sys.argv[i]:
        points = int(sys.argv[i].replace("--points=",""))
        args.append(sys.argv[i]); 
        i += 1
    elif sys.argv[i] in [ "--mass", "-m" ]:
        mass = float(sys.argv[i+1])
        args.append(sys.argv[i]); 
        args.append(sys.argv[i+1]); 
        i += 2
    elif sys.argv[i] in [ "-M", "--method" ]:
        method = sys.argv[i+1]
        args.append(sys.argv[i]); 
        args.append(sys.argv[i+1]); 
        i += 2
    elif sys.argv[i] == "--hadd":
        hadd = True
        i += 1
    else:
        args.append(sys.argv[i]); 
        i += 1
        
if points == 0: raise RuntimeError, "parallelScan requires that there be a --points=<n> or --points <n> option in the command line\n";
if jobs == 0:
    cpuinfo = open("/proc/cpuinfo","r")
    cores = sum([(1 if match("^processor\\b.*",l) else 0) for l in cpuinfo])
    if cores == 0: raise RuntimeError, "Cannot determine number of cores from /proc/cpuinfo, so I need a -j <n> option\n"; 
    if cores > 2 and match("(lxplus|cmslpc).*", os.environ['HOSTNAME']):
        jobs = cores/2
        print "Will run with %d jobs (half of the cores, for respect to other users)" % jobs
    else:
        jobs = cores
        print "Will run with %d jobs (one per core)" % jobs
if hadd:
    if not method: 
        print "Cannot understand what method of combine you are using, so cannot do hadd"
        exit()

workers = []
for j in xrange(jobs):
    start = int(ceil(points*j/float(jobs)))
    end   = int(ceil(points*(j+1)/float(jobs))-1)
    myargs = ["combine"] + args[:] + [ "-n", "%s.%d" % (name,j), "--firstPoint", str(start), "--lastPoint", str(end) ]
    print "spawning %s" % (" ".join(myargs))
    workers.append( subprocess.Popen(myargs) )

for w in workers:
    w.wait()

if hadd:
    if not method: 
        print "Cannot understand what method of combine you are using, so cannot do hadd"
    else:
        output = "higgsCombine%s.%s.mH%g.root" % (name,method,mass)
        input  = "  ".join(["higgsCombine%s.%d.%s.mH%g.root" % (name,j,method,mass) for j in xrange(jobs)])
        print "All workers done, doing the hadd to make %s out of higgsCombine%s.{0..%d}.%s.mH%g.root" % (output, name,jobs-1,method,mass)
        os.system("hadd -f %s %s" % (output,input))
else:
    print "All workers done, now you have to hadd the results yourself."
    
