#!/bin/bash
  
export R_HOME="/usr/local/lib/R"


gcc -Wall simu_discrete.c brent.c simulators.c condcp.c mle.c mleuser.c intervals.c -lRmath -lR -o main.out  

./main.out
