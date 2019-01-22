main_typeX.c: The main functions that print the data and MLEs.

r_object_typeX.c: Equivalent to main_typeX.c. But this function produce a R object.
compareMLE_typeX.r: Calling r_object_typeX.c from R

simulator_typeX.c: Functions that generate random samples.

derivative_typeX.c: This function is the equation that we want to find the root of.

brent.c: Downloaded online, to find out the root of the equation in derivative_typeX.c
brent.h: The corresponding 

overall_typeX.h: Contains the declarations of external variables and functions.