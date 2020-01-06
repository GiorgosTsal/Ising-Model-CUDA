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
	// I use two arrays, read from one and write to the other, then swap the pointers for the next iteration
	// Array to store local copy of G
	int *g_tmp = malloc(n*n * sizeof(int));

	// Temporary pointer to swap G and g_tmp
	int *swapg;

	// Variable to store the value of each moment
	double sum = 0;

	// The indices of the examined neighbors
	int idx_X, idx_Y;

	//Iterate k times
	for(int i = 0; i < k; i++)
	{
		//loop through every moment of G
		for(int x=0; x<n; x++)
			for(int y=0; y<n; y++)
			{
				sum = 0;

				//  Iterate through the moment's neighbors (k->X, l->Y axis)
				//	The neighborhood is a 5 × 5 window centered to each lattice point.
				for(int k=0; k<5; k++)
					for(int l=0; l<5; l++)
					{
						// Only edit the neighbors of the examined element
						if((k == 2) && (l == 2))
							continue;

						// Find the index of the examined neighbor 
						idx_X = (x + (k-2) + n) % n;
						idx_Y = (y + (l-2) + n) % n;

						// Calculate the new value
						sum += w[l*5 + k] * G[idx_Y*n + idx_X];
					}

				// If positive set state to 1
				// If negative set state to -1
				// If zero dont make changes
				if(sum > 0.001)
					g_tmp[y * n + x] = 1;
				else if(sum < -0.001)
					g_tmp[y * n + x] = -1;
				else
					g_tmp[y * n + x] = G[y * n + x];
			}

		// Swap pointers for next iteration
		swapg = G;
		G = g_tmp;
		g_tmp = swapg;
	}

	// At the last iteration, if the k is odd, G points to g_tmp and g_tmp points to G
	if(k%2 != 0)
		memcpy(g_tmp, G, n*n*sizeof(int));
}

int main()
{
	
	int n = 517;
	int k =11;

	// weight matrix 𝑤 
    double weights[] = {0.004, 0.016, 0.026, 0.016, 0.004,
                		0.016, 0.071, 0.117, 0.071, 0.016,
            			0.026, 0.117, 0    , 0.117, 0.026,
            			0.016, 0.071, 0.117, 0.071, 0.016,
            			0.004, 0.016, 0.026, 0.016, 0.004};

	// Open the binary file and write contents to G array
    FILE *fptr = fopen("conf/conf-init.bin","rb");
    if (fptr == NULL)
	{
        printf("Error opening file");
        exit(1);
    }
	// Array that will keep the init binary file info
	int *G = (int*)malloc(n*n * sizeof(int));
    fread(G, sizeof(int), n*n, fptr);
	fclose(fptr);


    // run ising procedure
    ising(G, weights, k, n);

	// name of conf file depending on k value(1,4,11)
	char filename[25];
	snprintf(filename, sizeof(filename), "conf/conf-%d.bin", k);

	// Compare returned data with the correct data 
	int *data = (int*)malloc(n*n * sizeof(int));
	int isWrong = 0;

	fptr = fopen(filename,"rb");
	fread(data, sizeof(int), n*n, fptr);
	fclose(fptr);
	for(int v = 0; v < n*n; v++)
		if(data[v] != G[v])
			isWrong = 1;


	if (!isWrong)
		printf("[k=%d] CORRECT\n", k);
	else
		printf("[k=%d] WRONG\n", k);

	// Free memory
    free(G);
	free(data);

    return 0;
}
