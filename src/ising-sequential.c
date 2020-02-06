/*
*       V0. Sequential: Simulation of an Ising model in two dimensions of size nxn for k iterations, staring from a uniform random initial state
*       Author:Tsalidis Georgios 27/12/2019
*       gtsalidis@ece.auth.gr
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/ising.h"

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
  int *g_tmp = new int[n*n]; //read from one

	// swap G and g_tmp pointeers for both matrixes
	int *swapg;

	// Variable to store the value of each moment
	double influence = 0;
 
  // Indexes of neibghors checked
	int idx_X, idx_Y;

	//flag to Terminate if no changesare made
	bool flag;
	//Iterate k times
	for(int i = 0; i < k; i++)
	{
		flag = false;
		//loop through G
		for(int x=0; x<n; x++)
			for(int y=0; y<n; y++)
			{
				influence = 0;

				// loop through the moment neighbors
				for(int X=0; X<5; X++)
					for(int Y=0; Y<5; Y++)
					{
						// find idx of checked point
						idx_X = (x + (X-2) + n) % n;
						idx_Y = (y + (Y-2) + n) % n;

						influence += *(w + Y*5 + X) * *(G +idx_Y*n + idx_X);
					}

		    //the value of the sign of influence If positive -> 1,If negative -> -1
				if(influence > 0.001){
					*(g_tmp + y*n + x) = 1;
					flag = true;
				}
				else if(influence < -0.001){
					*(g_tmp + y*n + x) = -1;
					flag = true;
				}
				else
       				  //remains the same
					*(g_tmp + y*n + x) = *(G + y*n + x);
				

		 }
		// Swap pointers for next iteration
		swapg = G;
		G = g_tmp;
		g_tmp = swapg;

		// Terminate model evolution if no changes were made
		if(!flag)
		{
			k = i+1;
			break;
		}

	
	}

  // Handle situation for odd 
	if(k%2 != 0)
		memcpy(g_tmp, G, n*n*sizeof(int));
 
}