#ifndef ISING_H
#define ISING_H



// Ising model evolution
/*!
  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]
  NOTE: Both matrices G and w are stored in row-major format.
*/
void ising( int *G, double *w, int k, int n);

#endif