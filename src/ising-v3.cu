/*
*       V3. GPU with multiple thread sharing common input moments
*       Author:Tsalidis Georgios 5/1/2020
*       gtsalidis@ece.auth.gr
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/ising.h"

#define CUDA_CHECK_ERROR() __cuda_check_errors(__FILE__, __LINE__)
#define LEN 2

#define BLOCK_X 128
#define BLOCK_Y 24
#define GRID_X 4
#define GRID_Y 4


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


__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int n, bool *flag);
bool evaluate(int *G1,int *G2, int n);


//kernel function used to calculate one thread with a block of moments
__global__ void ising_kernel(double* gpu_w, int* gpu_G, int* gpu_Gtmp, int n, bool *flag)
{
    int s_ncols = blockDim.x + 2*LEN;//num of  rows
 
    __shared__ double s_w[5*5];
    __shared__ int s_G[(BLOCK_X+2*LEN) * (BLOCK_Y+2*LEN)];


    for(int i=0; i<5*5; i++){
        *(s_w + i) = *(gpu_w+i);
    }

	// the value of each moment
	double influence;

	// moments x,y coordinates
	int x = blockDim.x*blockIdx.x + threadIdx.x;
    int y = blockDim.y*blockIdx.y + threadIdx.y;

	int s_x = threadIdx.x + LEN;
    int s_y = threadIdx.y + LEN;
 
    //thread increament
    int next_x = blockDim.x *gridDim.x ;
    int next_y = blockDim.y *gridDim.y ;

    //cordinates for neibghors on shared
    int n_x, n_y;

	// Indexes of neibghors checked
	int idx_x, idx_y;

	//each thread to compute a block of moments
	for(int i = x; i< n + LEN; i+= next_x)
	{
        for(int j = y; j< n + LEN; j+= next_y)
        {

            *(s_G + s_x*s_ncols + s_y) = *(gpu_G + ((i + n)%n)*n +  (j + n)%n);
         
         
			//For right and left
			if(threadIdx.x < LEN)
			{
				n_x = s_x;
				idx_x = (i + n)%n;

				for(int p=0; p<2; p++)
				{
					int count = (p-1)*LEN + p*blockDim.x;
					n_y = s_y + count;
					idx_y = (j + count + n)%n;
					s_G[n_x*s_ncols + n_y] = gpu_G[idx_x*n + idx_y];
				}
			}

			//For bot and top
			if(threadIdx.y < LEN)
			{
				n_y = s_y;
				idx_y = (j + n)%n;

				for(int p=0; p<2; p++)
				{
					int count = (p-1)*LEN + p*blockDim.y;
					n_x = s_x + count;
					idx_x = (i + count + n)%n;
					s_G[n_x*s_ncols + n_y] = gpu_G[idx_x*n + idx_y];
				}
			}

			//For corners
			if( (threadIdx.x < LEN) && (threadIdx.y<LEN) )
			{
				for(int p=0; p<4; p++)
				{
					int count_x = (p%2 - 1)*LEN + (p%2)*blockDim.y;
					n_x = s_x + count_x;
					idx_x = (i + count_x + n)%n;

					int count_y = ((p+3)%(p+1)/2 - 1)*LEN + ((p+3)%(p+1)/2)*blockDim.x;
					n_y = s_y + count_y;
					idx_y = (j + count_y + n)%n;

					s_G[n_x*s_ncols + n_y] = gpu_G[idx_x*n + idx_y];
				}
			}

			// Synchronize threads
			__syncthreads();

            if(i<n && j<n){
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
                     
                        influence += *(s_w + X*5 + Y) * *(s_G + (2+X+s_x)*s_ncols + (Y+s_x));
                    }

                //the value of the sign of influence If positive -> 1,If negative -> -1
                if(influence > 0.0001)
                {
                    *(gpu_Gtmp + i*n + j) = 1;
					*flag = true;
                }
                else if(influence < -0.0001)
                {
                    *(gpu_Gtmp + i*n + j) = -1;
					*flag = true;
                }
                else
                    //remains the same
                    *(gpu_Gtmp + i*n + j) = *(s_G + s_x*s_ncols + s_y);
            }
         __syncthreads();
        }
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

    //initialize blocks and threads => dim3 is an integer vector type based on uint3 
    //that is used to specify dimensions. When defining a variable of type dim3, any component left unspecified is initialized to 1

    dim3 block(BLOCK_X,BLOCK_Y); // blockDim
    dim3 grid(GRID_X,GRID_Y); // gridDim

	//flags for termination of no changes made
	bool flag;
	bool *gpu_flag;
	cudaMalloc(&gpu_flag, (size_t)sizeof(bool));

	//run for k iterations
	for(int i = 0; i < k; i++)
	{
		flag = false;
		cudaMemcpy(gpu_flag, &flag, (size_t)sizeof(bool), cudaMemcpyHostToDevice);

		//run kernel function to device
		ising_kernel<<< grid , block >>>(gpu_w, gpu_G, gpu_Gtmp, n, gpu_flag);
  
		//check for device errors
		CUDA_CHECK_ERROR ();
		
		//Synchronize 
		cudaDeviceSynchronize();

		//swap pointers 
		temp = gpu_G;
		gpu_G = gpu_Gtmp;
		gpu_Gtmp = temp;
		
		// Terminate model evolution if no changes were made
		cudaMemcpy(gpu_flag, &flag, (size_t)sizeof(bool), cudaMemcpyHostToDevice);
		if(flag)
		{
			break;
		}
	}

	cudaMemcpy(G, gpu_G, n*n*sizeof(int), cudaMemcpyDeviceToHost);

	// free GPU memory
	cudaFree(gpu_w);
	cudaFree(gpu_G);
	cudaFree(gpu_Gtmp);
}