#!/usr/bin/python

# Part b)
#


import sys
import time
import numpy as np

# beginning

if len(sys.argv) != 2:
  sys.exit("usage: " + sys.argv[0] + " <cov matrix file name>")

# fetch parameters
cov_filename = sys.argv[1]

print "reading data..."

# open files 

try:
    cov_file = open(cov_filename, 'r')
except IOError:
    print "Cannot open file %s\n" % cov_filename
    sys.exit("bye")


# import data from the returns file

N = int(cov_file.readline().split()[1])
q = np.ndarray(shape=(N, N))

cov_file.readline() # jump over the "matrix" word

for i in xrange(N):
    line = cov_file.readline().split()
    for j in xrange(N):
		q[i, j] = float(line[j]);

cov_file.close()

print "pca ..."
t = time.time()
eigval, eigvec = np.linalg.eig(q)
t = time.time() - t 
print "done"
print "PCA took " + str(round(1000.0*t)) + " ms"


sys.stdout.flush()

# end



