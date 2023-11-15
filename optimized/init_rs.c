/* Initialize a RS codec
 *
 * Copyright 2002 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */
#include "rs-common.h"

unsigned int nn=NN;

data_t ALPHA_TO[NN+1];
data_t INDEX_OF[NN+1];
data_t GENPOLY[E+1];

data_t recd[NN]; //Index_of data[]

/* Initialize a Reed-Solomon codec */
int init_rs(void)
{
  int i, j, sr,root;

  /* Generate Galois field lookup tables */
  INDEX_OF[0] = A0; /* log(zero) = -inf */
  ALPHA_TO[A0] = 0; /* alpha**-inf = 0 */

  sr = 1;
  for(i=0;i<NN;i++){
    INDEX_OF[sr] = i;
    ALPHA_TO[i] = sr;
    sr <<= 1;
    if(sr & (1<<MM))
      sr ^= GFPOLY;
    sr &= NN;
  }

  if(sr != 1)
  {
    /* field generator polynomial is not primitive! */
    return -1;
  }

  GENPOLY[0] = 1;
  for (i = 0,root=1; i < E; i++,root += 1) {
    GENPOLY[i+1] = 1;

    /* Multiply GENPOLY[] by  @**(root + x) */
    for (j = i; j > 0; j--){
      if (GENPOLY[j] != 0)
	GENPOLY[j] = GENPOLY[j-1] ^ ALPHA_TO[modnn(INDEX_OF[GENPOLY[j]] + root)];
      else
	GENPOLY[j] = GENPOLY[j-1];
    }
    /* GENPOLY[0] can never be zero */
    GENPOLY[0] = ALPHA_TO[modnn(INDEX_OF[GENPOLY[0]] + root)];
  }
  /* convert GENPOLY[] to index form for quicker encoding */
  for (i = 0; i <= E; i++)
    GENPOLY[i] = INDEX_OF[GENPOLY[i]];

  for(j=K+E;j<NN;j++)recd[j]=A0; //for decoder

  return 0;
}
