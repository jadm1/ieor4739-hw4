#!/bin/sh
cd runs
# ../bin/rpower ../data/q.txt -s 0 -q 1 -w 1 -r 10 -t 1e-7
../bin/covpca ../data/q.txt ../data/eigvals.txt ../data/eigvecs.txt -r 10 -t 1e-7