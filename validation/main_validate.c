#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/ising.h"

bool evaluate(int *G1,int *G2, int n ){
  for(int i=0;i<n*n;i++){
    if(G1[i]!=G2[i]){
      return false;
    }
  }
  return true;
}


int main() {
  
    int n = 517;
    //int k = 11;

    // weight matrix 𝑤 
    double weights[] = {0.004, 0.016, 0.026, 0.016, 0.004,
                		0.016, 0.071, 0.117, 0.071, 0.016,
            			0.026, 0.117, 0    , 0.117, 0.026,
            			0.016, 0.071, 0.117, 0.071, 0.016,
            			0.004, 0.016, 0.026, 0.016, 0.004};


    //Getting the initial situation of the lattice
    int *data =(int *)malloc((size_t)sizeof(int)*n*n);
    FILE *f = fopen("conf/conf-init.bin", "rb");
    fread(data,sizeof(int),n*n,f); // read bytes to our buffer

    int *stateInit =(int *)malloc((size_t)sizeof(int)*n*n);
    memcpy (stateInit, data, (size_t)sizeof(int)*n*n);

    //validate for k=1 
    memcpy (stateInit, data, (size_t)sizeof(int)*n*n);
    //Run ising proc
    ising(stateInit, weights,1, n);

    int *stateNxt_1 =(int *)malloc((size_t)sizeof(int)*n*n);
    FILE * fptr_1 = fopen("conf/conf-1.bin", "rb");
    fread(stateNxt_1,sizeof(int),n*n,fptr_1); // read file
    
    bool result_1 = evaluate(stateInit,stateNxt_1,n);
    if(result_1)
      printf("k=1: CORRECT \n");
    else
      printf("k=1: WRONG \n");
    
    free(stateNxt_1);

    //validate for k=4 
    memcpy (stateInit, data, (size_t)sizeof(int)*n*n);
    //Run ising proc
    ising(stateInit, weights,4, n);

    int *stateNxt_4 =(int *)malloc((size_t)sizeof(int)*n*n);
    FILE * fptr_4 = fopen("conf/conf-4.bin", "rb");
    fread(stateNxt_4,sizeof(int),n*n,fptr_4); // read file
   
    bool result_4=evaluate(stateInit,stateNxt_4,n);
    if(result_4)
      printf("k=4: CORRECT \n");
    else
      printf("k=4: WRONG \n");
  
    free(stateNxt_4);

    //validate for k=11 
    memcpy (stateInit, data, (size_t)sizeof(int)*n*n);
    //Run ising proc
    ising(stateInit, weights,11, n);

    int *stateNxt_11 =(int *)malloc((size_t)sizeof(int)*n*n);
    FILE * fptr_11 = fopen("conf/conf-11.bin", "rb");
    fread(stateNxt_11,sizeof(int),n*n,fptr_11); // read file
   
    bool result_11=evaluate(stateInit,stateNxt_11,n);
    if(result_11)
       printf("k=11: CORRECT \n");
    else
       printf("k=11: WRONG \n");

    free(stateNxt_11);

    return 0;
}