#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "LowLevel.h"
#include "timer.h"
#include "CCNT.h"

#include "rs-common.h"

static data_t ShuffleBuf[L];

static void ShuffleInit(void)
{
 for(data_t i=0;i<L;i++)ShuffleBuf[i]=i;
}

static void ShuffleProcess(void)
{
 for(data_t i=0;i<L;i++)
 {
  data_t p0=((data_t)rand())%L;
  data_t p1=((data_t)rand())%L;

  data_t t=ShuffleBuf[p0];
  ShuffleBuf[p0]=ShuffleBuf[p1];
  ShuffleBuf[p1]=t;
 }
}

static data_t block[L];
static data_t cblock[L];

static unsigned int t;

void main(void)
{
 LowLevel_Init();

 printf("RS GF(2^12)...\n");

 init_perfcounters(0,1);

 srand(get_cyclecount());

 ShuffleInit();

 if(init_rs())
 {
  printf("RS init fail!\n");;
  while(1);
 }

 Loop:

 memset(block,0,sizeof(block));

 /* Load block with random data */
 for(int i=0;i<K;i++)block[i]=((data_t)rand())&NN;

 /* Encode block */
 t=AVS_CNT0_REG;
 encode_rs(block);
 t=AVS_CNT0_REG-t;
 printf("Encode: %0.1lf FPS   ",(double)6000000.0/(double)t);

 /* Make copy, for compare */
 memcpy(cblock,block,sizeof(block));

 //Insert Errors
 ShuffleProcess();
 for(int i=0;i<(E/2)+0;i++)block[ShuffleBuf[i]]=(~block[ShuffleBuf[i]])&NN;

 /* Decode errored block */
 t=AVS_CNT0_REG;
 int r=decode_rs(block);
 t=AVS_CNT0_REG-t;

 if(r==-1)
 {
  printf("\n\nDecode Error!\n");
  while(1);
 }

 printf("Decode: %0.1lf FPS\n",(double)6000000.0/(double)t);

 if(memcmp(block,cblock,sizeof(block)))
 {
  printf("Decode Error!\n");
  while(1);
 }

 goto Loop;
}
