/*
*       V2. GPU with one thread computing a block of moments
*       Author:Tsalidis Georgios 5/1/2020
*       gtsalidis@ece.auth.gr
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/ising.h"

#define BLOCK_SIZE 128 // value usually chosen by tuning and hardware constraints
#define CUDA_CHECK_ERROR() __cuda_check_errors(__FILE__, __LINE__)

// See: http://codeyarns.com/2011/03/02/how-to-do-error-checking-in-cuda/
inline void
__cuda_check_errors (const char *filename, const int line_number)
{
  cudaError err = cudaDeviceSynchronize ();
  if (err != cudaSuccess)
    {
      printf ("CUDA error %i at %s:%i: %s\n",
          err, filename, line_number, cudaGetErrorString (err));
      exit (-1);
    }
}


__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int n);
bool evaluate(int *G1,int *G2, int n);


//kernel function used to calculate one thread with a block of moments
__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int n)
{
	//calculate thread_id
	int thread_id = blockIdx.x*blockDim.x + threadIdx.x;

	// the value of each moment
	double influence;

	// moments x,y coordinates
	int x, y;
	 
	// Indexes of neibghors checked
	int idx_x, idx_y;
	int next_thr = gridDim.x*blockDim.x;

	//each thread to compute a block of moments
	for(int thread = thread_id; thread<n*n; thread+= next_thr)
	{
		// moments x,y coordinates
		y = thread%n;
		x = thread/n;
		influence = 0;
		// loop through the moment neighbors
	    for(int X=0; X<5; X++)
	        for(int Y=0; Y<5; Y++)
	        {
				// skips the current iteration of the loop and continues with the next iteration.
	            if((X == 2) && (Y == 2))
	                continue;  
				
	            //find idx of checked point
	            idx_x = (x + (X-2) + n) % n;
	            idx_y = (y + (Y-2) + n) % n;

	            influence += *(gpu_w + X*5 + Y) * *(gpu_G +idx_x*n + idx_y);
	        }

	    //the value of the sign of influence If positive -> 1,If negative -> -1
		if(influence > 0.0001)
		{
			*(gpu_Gtmp + x*n + y) = 1;
		}
		else if(influence < -0.0001)
		{
			*(gpu_Gtmp + x*n + y) = -1;
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

	//dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
	//dim3 dimGrid(GRID_DIM_X,GRID_DIM_Y);
	int block = BLOCK_SIZE;
	int grid = (n + block -1)/block;

	//run for k iterations
	for(int i = 0; i < k; i++)
	{
		//run kernel function to device
		ising_kernel<<< grid , block >>>(gpu_w, gpu_G, gpu_Gtmp, n);
		
		//check for device errors
		CUDA_CHECK_ERROR ();
		
		//Synchronize 
		cudaDeviceSynchronize();

		//swap pointers 
		temp = gpu_G;
		gpu_G = gpu_Gtmp;
		gpu_Gtmp = temp;
	}

	cudaMemcpy(G, gpu_G, n*n*sizeof(int), cudaMemcpyDeviceToHost);

	// free GPU memory
	cudaFree(gpu_w);
	cudaFree(gpu_G);
	cudaFree(gpu_Gtmp);
}