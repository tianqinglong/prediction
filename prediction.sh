#!/bin/bash

export R_HOME=/usr/local/R/3.5.2/lib64/R
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/src/R-3.5.2/src/nmath/standalone:/usr/local/R/3.5.2/lib64/R/lib

gcc -Wall simu_continuous.c brent.c simulators.c condcp.c mle.c mleuser.c intervals.c -I /usr/local/R/3.5.2/lib64/R/include -L /usr/local/src/R-3.5.2/src/nmath/standalone -lRmath -lm -lR -L /usr/local/R/3.5.2/lib64/R/lib -lRblas -o main.out  

./main.out
