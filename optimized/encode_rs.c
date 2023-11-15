/* Reed-Solomon encoder
 * Copyright 2003, Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */
#include "rs-common.h"

void encode_rs(data_t *data)
{
 data_t feedback;

 data_t *parity=&data[K];

 CLEAR(parity,E);

 for(int i=0;i<K;i++)
 {
  feedback=data[i]^parity[0];

  if(feedback)
  {
   feedback=INDEX_OF[feedback];

   for(int j=0;j<E-1;j++)parity[j]  =parity[j+1]^ALPHA_TO[gf_add(GENPOLY[(E-1)-j],feedback)];
                         parity[E-1]=            ALPHA_TO[gf_add(GENPOLY[0]      ,feedback)];
  }
  else
  {
   for(int j=0;j<E-1;j++)parity[j]  =parity[j+1];
                         parity[E-1]=0;
  }

 }

}
