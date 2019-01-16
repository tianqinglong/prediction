#include <stdio.h>
#include <time.h>
#include <stdlib.h>
int main() {
    int i;
       srand(time(NULL)); /* set seed */
       for (i=1; i<=10; i++){
         printf("%f ", rand() /  (double) RAND_MAX );
}
return 0; }
