@echo off
cd /d %~dp0

echo --- Myopt ---
bin\myopt data\myoinput.txt data\x.txt -m 1000 -t 1e-6
echo ---------------
echo --- P&L simulations ---
bin\pfsimul data\x.txt data\p.txt -ov runs\pf_values.txt -or runs\pf_returns.txt -q 1000000 -w 4 -B 1000000000 -p 504 -v

pause

