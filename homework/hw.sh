#!/bin/sh
cd runs
../bin/rpower ../data/russell_1000_cov.txt -s 0 -q 1 -w 1 -r 10 -t 1e-6 # only 1 unperturbed job
#../bin/rpower ../data/russell_1000_cov.txt -s 10 -q 8 -w 4 -r 10 -t 1e-6
