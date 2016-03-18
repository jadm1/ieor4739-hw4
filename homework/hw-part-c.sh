#!/bin/sh
cd runs
../bin/myopt ../data/myoinput.txt ../data/x.txt -m 100 -t 1e-6 # get an optimal portfolio allocation
time ../bin/pfsimul ../data/x.txt ../data/p.txt pf_values.txt pf_returns.txt -q 100000 -w 4

