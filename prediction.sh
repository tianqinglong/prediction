#!/bin/bash

export R_HOME="/usr/local/Cellar/r/3.5.2_2/lib/R"

gcc -Wall simu_continuous.c brent.c simulators.c condcp.c mle.c mleuser.c intervals.c -lRmath -lm -lR -o main.out  

./main.out