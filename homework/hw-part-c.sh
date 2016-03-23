#!/bin/sh

echo --- Myopt ---
bin/myopt data/myoinput.txt data/x.txt -m 100 -t 1e-6 # get an optimal portfolio allocation
echo ---------------
echo --- PnL simulations ---
bin/pfsimul data/x.txt data/p.txt pf_values.txt pf_returns.txt -v -q 10000 -w 4
