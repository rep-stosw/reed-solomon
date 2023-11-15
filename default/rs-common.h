/* Stuff common to all the general-purpose Reed-Solomon codecs
 * Copyright 2004 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */

#ifndef _RS_COMMON_H_
#define _RS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

//#define UNNORMALIZED

#undef MIN
#define	MIN(a,b)	((a) < (b) ? (a) : (b))

#define MODNN(x) modnn(x)

#define MM 12
#define NN ((1<<MM)-1)
#define A0  NN

#define K   3189
#define E   779
#define L   (K+E)
#define PAD (NN-L)

#define NROOTS E

#define GFPOLY 016533

typedef unsigned short data_t;

extern data_t ALPHA_TO[NN+1];
extern data_t INDEX_OF[NN+1];
extern data_t GENPOLY[NROOTS+1];

static inline int modnn(int x)
{
  while (x >= NN) {
    x -= NN;
    x = (x >> MM) + (x & NN);
  }
  return x;
}

int init_rs(void);
void encode_rs(data_t *data);
int decode_rs(data_t *data);

#ifdef __cplusplus
}
#endif

#endif
