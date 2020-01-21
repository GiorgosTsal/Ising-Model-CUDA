/*
*       V0. Sequential 27/12/2019
*       Author:Tsalidis Georgios
*       gtsalidis@ece.auth.gr
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Ising model evolution
/*
  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]
  NOTE: Both matrices G and w are stored in row-major format.
*/
void ising(int *G, double *w, int k, int n)
{
	// Use two arrays, read from one and write to the other, then swap the pointers for the next iteration
  int *g_tmp = new int[n*n];

	// swap G and g_tmp pointeers for both matrixes
	int *swapg;

	// Variable to store the value of each moment
	double total = 0;

	// The indices of the checked neighbors
	int idx_X, idx_Y;

	//Iterate k times
	for(int i = 0; i < k; i++)
	{
		//scan the lattice-loop through G
		for(int x=0; x<n; x++)
			for(int y=0; y<n; y++)
			{
				total = 0;

				//  Iterate through the moments neighbors (k->X, l->Y axis) to compute the next sping
				//	The neighborhood is a 5 × 5 window centered to each lattice point.
				for(int k=0; k<5; k++)
					for(int l=0; l<5; l++)
					{
						// Only edit the neighbors of the examined element
						if((k == 2) && (l == 2))
							continue;

						// Find the index of the examined neighbor with modulus
						idx_X = (x + (k-2) + n) % n;
						idx_Y = (y + (l-2) + n) % n;

						// Calculate the new influence value
						total += w[l*5 + k] * G[idx_Y*n + idx_X];
					}

				// If positive set state to 1
				// If negative set state to -1
				// If zero dont make changes
				if(total > 0.001)
					g_tmp[y * n + x] = 1;
				else if(total < -0.001)
					g_tmp[y * n + x] = -1;
				else
					g_tmp[y * n + x] = G[y * n + x];
			}

		// Swap pointers for next iteration
		swapg = G;
		G = g_tmp;
		g_tmp = swapg;
	}

  // Handle situation if k is odd at the last iteration
	if(k%2 != 0)
		memcpy(g_tmp, G, n*n*sizeof(int));
}

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
    int k = 11;

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