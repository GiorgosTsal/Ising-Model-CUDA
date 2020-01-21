
/*
*       V1. GPU with one thread per moment 
*       Author:Tsalidis Georgios 2/1/2020
*       gtsalidis@ece.auth.gr
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK_SIZE 128 //threads per block (^2)

__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int* c, int n);
bool evaluate(int *G1,int *G2, int n);


//kernel function used to calculate one moment per thread
__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int* c, int n)
{
	//calculate thread_id
	int thread_id = blockIdx.x*blockDim.x + threadIdx.x;

	// the value of each moment
	double influence;

	// moments x,y coordinates
	int y = thread_id%n;
	int x = thread_id/n;
	
	// Indexes of neibghors checked
	int idx_x, idx_y;

	if( thread_id < n*n )
	{
		// loop through the moment neighbors
	    for(int X=0; X<5; X++)
	        for(int Y=0; Y<5; Y++)
	        {
	            //find idx of checked point
	            idx_x = (x + (X-2) + n) % n;
	            idx_y = (y + (Y-2) + n) % n;

	            influence += *(gpu_w + X*5 + Y) * *(gpu_G +idx_x*n + idx_y);
	        }

	    //the value of the sign of influence If positive -> 1,If negative -> -1
		if(influence > 0.0001)
		{
			*(gpu_Gtmp + x*n + y) = 1;
			*c = 1;
		}
		else if(influence < -0.0001)
		{
			*(gpu_Gtmp + x*n + y) = -1;
			*c = 1;
		}
	    else
					//remains the same
	        *(gpu_Gtmp + x*n + y) = *(gpu_G + x*n + y);
	}
}


void ising(int *G, double *w, int k, int n)
{
   
	double *gpu_w;
	int *gpu_G;
	
	// allocate weight array and G array
	cudaMalloc(&gpu_w, 5*5*sizeof(double));
	cudaMalloc(&gpu_G, n*n*sizeof(int));
	
	
	//transfer data to device(GPU)
	cudaMemcpy(gpu_w, w, 5*5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_G, G, n*n*sizeof(int), cudaMemcpyHostToDevice);

	//GPU array to store the updated values
	int *gpu_Gtmp;
	cudaMalloc(&gpu_Gtmp, n*n*sizeof(int));

	// gpu_G with gpu_Gtmp pointer swap
	int *temp;

	int threads = BLOCK_SIZE;
	int blocks = (n*n + threads - 1)/threads;

	//vars to store changes 
	int c;
	int *gpu_c;

	cudaMalloc(&gpu_c, (size_t)sizeof(int));

	//run for k iterations
	for(int i = 0; i < k; i++)
	{
		c = 0;

		cudaMemcpy(gpu_c, &c, (size_t)sizeof(int), cudaMemcpyHostToDevice);
		
		//run kernel function to device
		ising_kernel<<< blocks , threads >>>(gpu_w, gpu_G, gpu_Gtmp, gpu_c, n);

		//Synchronize 
		cudaDeviceSynchronize();

		//swap pointers 
		temp = gpu_G;
		gpu_G = gpu_Gtmp;
		gpu_Gtmp = temp;

		cudaMemcpy(&c, gpu_c,  (size_t)sizeof(int), cudaMemcpyDeviceToHost);
		
		//if no changes made, break
		if(c == 0)
			break;
	}

	cudaMemcpy(G, gpu_G, n*n*sizeof(int), cudaMemcpyDeviceToHost);

	// free GPU memory
	cudaFree(gpu_w);
	cudaFree(gpu_G);
	cudaFree(gpu_Gtmp);
}

//function for validation 
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