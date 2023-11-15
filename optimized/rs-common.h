/* Stuff common to all the general-purpose Reed-Solomon codecs
 * Copyright 2004 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */
#ifndef _RS_COMMON_H_
#define _RS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#define	MIN(a,b) ((a)<(b)?(a):(b))

#define MM 12
#define NN ((1<<MM)-1)
#define A0  NN

#define K   3189
#define E   779
#define L   (K+E)
#define PAD (NN-L)

#if PAD < 0
#error "Invalid Parameters: K,E"
#endif

/*
for GF(2^12), octal

       010123,  015647,  016533,  016047,  011015,  014127,
       017673,  013565,  015341,  015053,  015621,  015321,
       011417,  013505
*/
#define GFPOLY 016533

typedef unsigned short data_t;

extern data_t ALPHA_TO[NN+1];
extern data_t INDEX_OF[NN+1];
extern data_t GENPOLY[E+1];

extern data_t recd[NN]; //for decoder

extern unsigned int nn;

static inline unsigned int modnn(unsigned int x)
{
 while(x>=nn)
 {
  x-=nn;
  x=(x>>MM)+(x&nn);
 }
 return x;
}

static inline unsigned int modnnF(unsigned int x)
{
 if(x<NN)return x;
 asm volatile ("" ::: "memory");
 return x-NN;
}

static inline unsigned int gf_sub(unsigned int x,unsigned int y)
{
 y=x-y;
 return y+(y>x)*NN;
}

#define gf_add(x,y) gf_sub(x,NN-(y))

#if 0

#include <string.h>

#define	CLEAR(a,n)      memset(a,0,(n)*sizeof((a)[0]))
#define COPY(a,b,n)     memcpy(a,b,(n)*sizeof((a)[0]))
#define	COPYDOWN(a,b,n) memmove(a,b,(n)*sizeof((a)[0]))

#else

#define	CLEAR(a,n) {\
	int ci;\
	for(ci=(n)-1;ci >=0;ci--)\
		(a)[ci] = 0;\
	}

#define	COPY(a,b,n) {\
	int ci;\
	for(ci=(n)-1;ci >=0;ci--)\
		(a)[ci] = (b)[ci];\
	}
#define	COPYDOWN(a,b,n) {\
	int ci;\
	for(ci=(n)-1;ci >=0;ci--)\
		(a)[ci] = (b)[ci];\
	}

#endif

int init_rs(void);
void encode_rs(data_t *data);
int decode_rs(data_t *data);

#ifdef __cplusplus
}
#endif

#endif
