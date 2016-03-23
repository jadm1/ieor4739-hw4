@echo off
cd /d %~dp0
python python\compute_vfqd.py data\returns.txt data\family.txt data\v.txt data\f.txt data\q.txt data\d.txt data\myoinput.txt
pause
