#!/bin/sh
cd runs
../bin/myopt ../data/myoinput.txt ../data/x.txt -m 100 -t 1e-7 # get an optimal portfolio allocation
../bin/pfsimul

