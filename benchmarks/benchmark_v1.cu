/*
*       Benchmarking GPU with one thread per moment 
*       Author:Tsalidis Georgios 2/1/2020
*       gtsalidis@ece.auth.gr
*/

#include  <time.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/ising.h"

struct timespec start, finish;
double elapsed;


static int n[]={50,300,600,900,1400,1800,2500,3800,6000};//9
static int k[]={1, 5, 10 , 30 ,60, 100};//6


int main() {

  printf("=====================BENCHMARKING CUDA V1 START=====================\n");

  double weights[]={ 0.004,  0.016,  0.026,  0.016,   0.004,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.026,  0.117,    0  ,  0.117,   0.026,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.004,  0.016,  0.026,  0.016,   0.004};


  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 6; j++) {
      FILE *pointerToFile;
      int * sample;

      sample=(int *)malloc((size_t)n[i]*n[i]*sizeof(int));
      for (size_t ii = 0; ii < n[i]*n[i]; ii++) {
        sample[ii]=(rand()%2 );
        if(sample[ii]==0){
          sample[ii]=-1;
        }
      }


      clock_gettime(CLOCK_MONOTONIC, &start);
      ising(sample, weights,k[j], n[i]);
      clock_gettime(CLOCK_MONOTONIC, &finish);
      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

      pointerToFile=fopen("benchmarks/results_v1.csv","a");
      fprintf(pointerToFile,"%d,%d,%lf\n",n[i],k[j],elapsed);
      printf("Ising model evolution for n=%d, k=%d ,took %lf seconds! \n",n[i],k[j], elapsed );
      free(sample);
    }

  }
  printf("\n");

  printf("=====================BENCHMARKING CUDA V1 END=====================\n");


}