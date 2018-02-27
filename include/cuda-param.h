#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>

#include "param.h"

#define RANDOM_SEED 1234
//#define RANDOM_MAX  10000

//#define MAX_TREADS_PER_BLOCK (int)1024
//#define MULTIPROCESSOR (int)15
//#define NUM_LOADS (1+(N/(int)(MULTIPROCESSOR*MAX_TREADS_PER_BLOCK)))


//#define BLOCKS ((int)NUM_LOADS*(int)MULTIPROCESSOR)
//#define THREADS ((N/(int)(BLOCKS))+1)

#define BLOCKS 55
#define THREADS 728


