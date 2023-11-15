/* General purpose Reed-Solomon decoder
 * Copyright 2003 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */
#include <stdio.h>

#include "rs-common.h"

#include "timer.h"

//#define SYNDROME_CLASSIC

#ifndef SYNDROME_CLASSIC

#define N1 13  //4-stage decomposition
#define N2  9
#define N3  7
#define N4  5

#define A 2835 //4D interleaver coeffs. (Chinese Remainder Theorem solution)
#define B  910
#define C 1170
#define D 3276

#define GAMMA(i1,i2,i3,i4) ((((i1)*A)+((i2)*B)+((i3)*C)+((i4)*D))%NN)             /* Mapping indexes i1,i2,i3,i4 to common index i */

#define CONDITION(x1,x2,x3,x4) if(GAMMA(x1,x2,x3,x4)&&(GAMMA(x1,x2,x3,x4)<(E+1))) /* We need only syndromes from 1 to E */

#define GAMMA_1(i1) ((i1)*A)
#define GAMMA_2(i2) ((i2)*B)
#define GAMMA_3(i3) ((i3)*C)
#define GAMMA_4(i4) ((i4)*D)

static unsigned int i1,i2,i3,i4;
static unsigned int j1,j2,j3,j4;

static data_t s1[NN+1];
static data_t s2[NN+1];
static data_t s3[NN+1];

#endif

int decode_rs(data_t *data)
{
  static data_t lambda[E+1];
  static data_t s[E+1];	     /* Err Locator poly and syndrome poly */
  static data_t b[E+1];
  static data_t t[E+1];
  static data_t omega[E+1];
  static data_t root[E];
  static data_t reg[E+1];
  static data_t Q[K+E];

  int deg_lambda, el, deg_omega;
  int i, j, r,k;
  int syn_error, count;

  data_t u,q,tmp,num1,den,discr_r;

  for(j=0;j<K+E;j++)recd[j]=INDEX_OF[data[(K+E-1)-j]]; //re-ordered & indexed data for optimization

//unsigned int tt=AVS_CNT0_REG;

//Classic ---------------------------------------------------------------------------------

#ifdef SYNDROME_CLASSIC

  syn_error=0;

  for(i=1;i<=E;i++)
  {
   tmp=0;

   for(j=0;j<K+E;j++)if(recd[j]!=A0)tmp^=ALPHA_TO[modnn(recd[j]+(i*j))];

   syn_error|=tmp;
   s[i]=tmp;
  }

#else

//Multi-Stage -----------------------------------------------------------------------------

//Stage 1
 for(i1=0;i1<N1;i1++)for(i2=0;i2<N2;i2++)for(i3=0;i3<N3;i3++)for(j4=0;j4<N4;j4++)
 {
  tmp=0;
  for(i4=0;i4<N4;i4++)if(recd[GAMMA(i1,i2,i3,i4)]!=A0)tmp^=ALPHA_TO[modnn(recd[GAMMA(i1,i2,i3,i4)]+GAMMA_4((i4*j4)%N4))];
  s1[GAMMA(i1,i2,i3,j4)]=tmp;
 }

 for(i=0;i<NN;i++)s1[i]=INDEX_OF[s1[i]];

//Stage 2
 for(i1=0;i1<N1;i1++)for(i2=0;i2<N2;i2++)for(j3=0;j3<N3;j3++)for(j4=0;j4<N4;j4++)
 {
  tmp=0;
  for(i3=0;i3<N3;i3++)if(s1[GAMMA(i1,i2,i3,j4)]!=A0)tmp^=ALPHA_TO[modnn(s1[GAMMA(i1,i2,i3,j4)]+GAMMA_3((i3*j3)%N3))];
  s2[GAMMA(i1,i2,j3,j4)]=tmp;
 }

 for(i=0;i<NN;i++)s2[i]=INDEX_OF[s2[i]];

//Stage 3
 for(i1=0;i1<N1;i1++)for(j2=0;j2<N2;j2++)for(j3=0;j3<N3;j3++)for(j4=0;j4<N4;j4++)
 {
  tmp=0;
  for(i2=0;i2<N2;i2++)if(s2[GAMMA(i1,i2,j3,j4)]!=A0)tmp^=ALPHA_TO[modnn(s2[GAMMA(i1,i2,j3,j4)]+GAMMA_2((i2*j2)%N2))];
  s3[GAMMA(i1,j2,j3,j4)]=tmp;
 }

 for(i=0;i<NN;i++)s3[i]=INDEX_OF[s3[i]];

 syn_error=0;

//Stage 4
 for(j4=0;j4<N4;j4++)for(j3=0;j3<N3;j3++)for(j2=0;j2<N2;j2++)for(j1=0;j1<N1;j1++)CONDITION(j1,j2,j3,j4)
 {
  tmp=0;
  for(i1=0;i1<N1;i1++)if(s3[GAMMA(i1,j2,j3,j4)]!=A0)tmp^=ALPHA_TO[modnn(s3[GAMMA(i1,j2,j3,j4)]+GAMMA_1((i1*j1)%N1))];
  s[GAMMA(j1,j2,j3,j4)]=tmp;

  syn_error|=tmp; /* set flag if non-zero syndrome => error */
 }

//------------------------------------------------------------------------------------------------------

#endif

// tt=AVS_CNT0_REG-tt;
// printf("Syndromes: %0.1lf FPS   ",(double)6000000.0/(double)tt);


  if (!syn_error)
  {
   /* if syndrome is zero, data[] is a codeword and there are no errors to correct. So return data[] unmodified */
   return 0;
  }

  /* Convert syndromes to index form, checking for nonzero condition */
  for(i=1;i<=E;i++)s[i]=INDEX_OF[s[i]];

  CLEAR(&lambda[1],E);
  lambda[0] = 1;

  for(i=1;i<E+1;i++)b[i]=A0;
  b[0]=0;

  /* Begin Berlekamp-Massey algorithm to determine error locator polynomial */
  r=0;
  el=0;
  while(++r<=E) /* r is the step number */
  {
   /* Compute discrepancy at the r-th step in poly-form */
   discr_r=0;

   for(i=0;i<r;i++)if(lambda[i]&&(s[r-i]!=A0))discr_r^=ALPHA_TO[gf_add(INDEX_OF[lambda[i]],s[r-i])];

   discr_r=INDEX_OF[discr_r]; /* Index form */

   if(discr_r==A0)
   {
    /* 2 lines below: B(x) <-- x*B(x) */
    COPYDOWN(&b[1],b,E);
    b[0]=A0;
   }
   else
   {
    /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
    t[0]=lambda[0];

    for(i=0;i<E;i++)
    {
     if(b[i]!=A0)t[i+1]=lambda[i+1]^ALPHA_TO[gf_add(discr_r,b[i])];
     else        t[i+1]=lambda[i+1];
    }

    if(el*2<r)
    {
     el=r-el;

     /* 2 lines below: B(x) <-- inv(discr_r) lambda(x) */
     for(i=0;i<=E;i++)b[i]=(!lambda[i])?A0:gf_sub(INDEX_OF[lambda[i]],discr_r);
    }
    else
    {
     /* 2 lines below: B(x) <-- x*B(x) */
     COPYDOWN(&b[1],b,E);
     b[0]=A0;
    }

    COPY(lambda,t,(E+1));
   }
  } //while

  /* Convert lambda to index form and compute deg(lambda(x)) */
  deg_lambda=0;

  for(i=0;i<E+1;i++)
  {
   lambda[i]=INDEX_OF[lambda[i]];
   if(lambda[i]!=A0)deg_lambda=i;
  }

  /* Find roots of the error locator polynomial by Chien search */
  COPY(&reg[1],&lambda[1],E);
  count=0;                    /* Number of roots of lambda(x) */

//Chien search Optimized ------------------------------------------------

//unsigned int tt=AVS_CNT0_REG;

  for(i=0;i<K+E;i++)Q[i]=1;

  if(PAD)for(j=deg_lambda;j>0;j--)if(reg[j]!=A0)reg[j]=modnn(reg[j]+(j*PAD));

  for(j=deg_lambda;j>0;j--)if(reg[j]!=A0)for(i=K+E-1;i>=0;i--)
  {
   reg[j]=modnnF(reg[j]+j);
   Q[i]^=ALPHA_TO[reg[j]];
  }

  for(i=K+E-1;i>=0;i--)if(!Q[i])
  {
   root[count]=NN-i;             /* store root (index-form) */
   if(++count==deg_lambda)break; /* If we've already found max possible roots, abort the search to save time */
  }

// tt=AVS_CNT0_REG-tt;
// printf("Chien: %0.1lf FPS   ",(double)6000000.0/(double)tt);

//-----------------------------------------------------------------------

  if(deg_lambda!=count)
  {
   /* deg(lambda) unequal to number of roots => uncorrectable error detected */
   return -1;
  }

  /*
   * Compute err evaluator poly omega(x) = s(x)*lambda(x) (modulo
   * x**E). in index form. Also find deg(omega).
   */
  deg_omega=deg_lambda-1;
  for(i=0;i<=deg_omega;i++)
  {
   tmp=0;
   for(j=i;j>=0;j--)if((s[i-j+1]!=A0)&&(lambda[j]!=A0))tmp^=ALPHA_TO[gf_add(s[i-j+1],lambda[j])];
   omega[i]=INDEX_OF[tmp];
  }

  /*
   * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
   * inv(X(l))**(FCR-1) and den = lambda_pr(inv(X(l))) all in poly-form
   */
  for (j = count-1; j >=0; j--)
  {
   num1=0;
   for(i=deg_omega;i>=0;i--)if(omega[i]!=A0)num1^=ALPHA_TO[modnn(omega[i]+i*root[j])];

   den=0;
   /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
   for(i=MIN(deg_lambda,E-1)&~1;i>=0;i-=2)if(lambda[i+1]!=A0)den^=ALPHA_TO[modnn(lambda[i+1]+i*root[j])];

   /* Apply error to data */
   if(num1&&(root[j]>PAD))data[root[j]-(PAD+1)]^=ALPHA_TO[gf_sub(INDEX_OF[num1],INDEX_OF[den])];
  }

 return count;
}
