/* General purpose Reed-Solomon decoder
 * Copyright 2003 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */

/* The guts of the Reed-Solomon decoder, meant to be #included
 * into a function body with the following typedefs, macros and variables supplied
 * according to the code parameters:

 * data_t - a typedef for the data symbol
 * data_t data[] - array of NN data and parity symbols to be corrected in place
 * retval - an integer lvalue into which the decoder's return code is written
 * NROOTS - the number of roots in the RS code generator polynomial,
 *          which is the same as the number of parity symbols in a block.
            Integer variable or literal.
 * NN - the total number of symbols in a RS block. Integer variable or literal.
 * PAD - the number of pad symbols in a block. Integer variable or literal.
 * ALPHA_TO - The address of an array of NN elements to convert Galois field
 *            elements in index (log) form to polynomial form. Read only.
 * INDEX_OF - The address of an array of NN elements to convert Galois field
 *            elements in polynomial form to index (log) form. Read only.
 * MODNN - a function to reduce its argument modulo NN. May be inline or a macro.

 * The memset(), memmove(), and memcpy() functions are used. The appropriate header
 * file declaring these functions (usually <string.h>) must be included by the calling
 * program.
 */
#include <string.h>
#include <stdio.h>

#include "rs-common.h"

#include "timer.h"

int decode_rs(data_t *data)
{
  static data_t lambda[NROOTS+1];
  static data_t s[NROOTS];	/* Err+Eras Locator poly
				 * and syndrome poly */
  static data_t b[NROOTS+1];
  static data_t t[NROOTS+1];
  static data_t omega[NROOTS+1];
  static data_t root[NROOTS];
  static data_t reg[NROOTS+1];
  static data_t loc[NROOTS];

  int deg_lambda, el, deg_omega;
  int i, j, r,k;
  int syn_error, count;

  data_t u,q,tmp,num1,den,discr_r;

  int retval;

  /* form the syndromes; i.e., evaluate data(x) at roots of g(x) */
  for(i=0;i<NROOTS;i++)
    s[i] = data[0];

//unsigned int tt=AVS_CNT0_REG;

  for(j=1;j<NN-PAD;j++)
  {
   for(i=0;i<NROOTS;i++)
   {

    if(s[i] == 0)s[i] = data[j];
    else         s[i] = data[j] ^ ALPHA_TO[MODNN(INDEX_OF[s[i]] + i + 1  )];

   }
  }

// tt=AVS_CNT0_REG-tt;
// printf("Syndromes: %0.1lf FPS   ",(double)6000000.0/(double)tt);


  /* Convert syndromes to index form, checking for nonzero condition */
  syn_error = 0;
  for(i=0;i<NROOTS;i++){
    syn_error |= s[i];
    s[i] = INDEX_OF[s[i]];
  }

  if (!syn_error) {
    /* if syndrome is zero, data[] is a codeword and there are no
     * errors to correct. So return data[] unmodified
     */
    count = 0;
    goto finish;
  }
  memset(&lambda[1],0,NROOTS*sizeof(lambda[0]));
  lambda[0] = 1;

  for(i=0;i<NROOTS+1;i++)b[i] = INDEX_OF[lambda[i]];

  /*
   * Begin Berlekamp-Massey algorithm to determine error+erasure
   * locator polynomial
   */
  r = 0; //no_eras;
  el = 0; //no_eras;
  while (++r <= NROOTS) {	/* r is the step number */
    /* Compute discrepancy at the r-th step in poly-form */
    discr_r = 0;
    for (i = 0; i < r; i++){
      if ((lambda[i] != 0) && (s[r-i-1] != A0)) {
	discr_r ^= ALPHA_TO[MODNN(INDEX_OF[lambda[i]] + s[r-i-1])];
      }
    }
    discr_r = INDEX_OF[discr_r];	/* Index form */
    if (discr_r == A0) {
      /* 2 lines below: B(x) <-- x*B(x) */
      memmove(&b[1],b,NROOTS*sizeof(b[0]));
      b[0] = A0;
    } else {
      /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
      t[0] = lambda[0];
      for (i = 0 ; i < NROOTS; i++) {
	if(b[i] != A0)
	  t[i+1] = lambda[i+1] ^ ALPHA_TO[MODNN(discr_r + b[i])];
	else
	  t[i+1] = lambda[i+1];
      }
      if (2 * el <= r - 1) {
	el = r - el;
	/*
	 * 2 lines below: B(x) <-- inv(discr_r) *
	 * lambda(x)
	 */
	for (i = 0; i <= NROOTS; i++)
	  b[i] = (lambda[i] == 0) ? A0 : MODNN(INDEX_OF[lambda[i]] - discr_r + NN);
      } else {
	/* 2 lines below: B(x) <-- x*B(x) */
	memmove(&b[1],b,NROOTS*sizeof(b[0]));
	b[0] = A0;
      }
      memcpy(lambda,t,(NROOTS+1)*sizeof(t[0]));
    }
  }

  /* Convert lambda to index form and compute deg(lambda(x)) */
  deg_lambda = 0;
  for(i=0;i<NROOTS+1;i++){
    lambda[i] = INDEX_OF[lambda[i]];
    if(lambda[i] != A0)
      deg_lambda = i;
  }
  /* Find roots of the error+erasure locator polynomial by Chien search */
  memcpy(&reg[1],&lambda[1],NROOTS*sizeof(reg[0]));
  count = 0;		/* Number of roots of lambda(x) */

//unsigned int tt=AVS_CNT0_REG;

  for(i=NN-1;i>=0;i--)
  {
   q = 1; /* lambda[0] is always 0 */

   for(j=deg_lambda;j>0;j--)
   {
    if(reg[j]!=A0)
    {
     reg[j]=MODNN(reg[j]+j);
     q^=ALPHA_TO[reg[j]];
    }
   }

   if(q)continue; /* Not a root */

   /* store root (index-form) and error location number */
   root[count]=NN-i;
   loc[count]= NN-i-1;

   /* If we've already found max possible roots,
    * abort the search to save time
   */
   if(++count==deg_lambda)break;
  }

// tt=AVS_CNT0_REG-tt;
// printf("Chien: %0.1lf FPS   ",(double)6000000.0/(double)tt);

  if (deg_lambda != count) {
    /*
     * deg(lambda) unequal to number of roots => uncorrectable
     * error detected
     */
    count = -1;
    goto finish;
  }

  /*
   * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
   * x**NROOTS). in index form. Also find deg(omega).
   */
  deg_omega = deg_lambda-1;
  for (i = 0; i <= deg_omega;i++){
    tmp = 0;
    for(j=i;j >= 0; j--){
      if ((s[i - j] != A0) && (lambda[j] != A0))
	tmp ^= ALPHA_TO[MODNN(s[i - j] + lambda[j])];
    }
    omega[i] = INDEX_OF[tmp];
  }

  /*
   * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
   * inv(X(l))**(FCR-1) and den = lambda_pr(inv(X(l))) all in poly-form
   */
  for (j = count-1; j >=0; j--) {
    num1 = 0;
    for (i = deg_omega; i >= 0; i--) {
      if (omega[i] != A0)
	num1  ^= ALPHA_TO[MODNN(omega[i] + i * root[j])];
    }

    den = 0;
    
    /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
    for (i = MIN(deg_lambda,NROOTS-1) & ~1; i >= 0; i -=2) {
      if(lambda[i+1] != A0)
	den ^= ALPHA_TO[MODNN(lambda[i+1] + i * root[j])];
    }

    /* Apply error to data */
    if (num1 != 0 && loc[j] >= PAD) {
      data[loc[j]-PAD] ^= ALPHA_TO[MODNN(INDEX_OF[num1] + NN - INDEX_OF[den])];
    }
  }

 finish:
  retval = count;
  
  return retval;
}
