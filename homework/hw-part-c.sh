#!/bin/sh

echo --- Myopt ---
bin/myopt data/myoinput.txt data/x.txt -m 1000 -t 1e-6 # get an optimal portfolio allocation
echo ---------------
echo --- PnL simulations ---
bin/pfsimul data/x.txt data/p.txt -ov runs/pf_values.txt -or runs/pf_returns.txt -q 1000000 -w 4 -b 1000000000 -p 504 -v
