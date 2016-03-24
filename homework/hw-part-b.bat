@echo off
cd /d %~dp0

echo ---Power method PCA using C---
echo for 10 eigen values
bin\covpca data\q.txt data\eigvals.txt data\eigvecs.txt -r 10 -t 1e-6 | grep PCA
echo ---------------
echo ---Numpy PCA---
echo for 947 eigen values
python python\covpca.py data\q.txt

