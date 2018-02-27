#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <curand.h>
#include <curand_kernel.h>

#include "param.h"
#include "particle.h"
#include "cuda-param.h"

#define TWO_PI 6.283185f

__global__ void G_initialize_random_generator(unsigned int seed, curandState_t *state)
{
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  curand_init(seed, idx, 0, &state[idx]);
};

__global__ void G_initialize(particle *G_boid,  curandState_t *G_random_state, unsigned int *G_random_number)
{
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  float angle;
  float auxiliar_rand;
  //inicializacao das posicoes e velocidades das particulas
  if(idx<N1)
    {
      angle = curand_uniform(&G_random_state[idx]);
      angle = TWO_PI * angle;
      G_boid[idx].vx = V1 * cos(angle);
      G_boid[idx].vy = V1 * sin(angle);
      auxiliar_rand = curand_uniform(&G_random_state[idx]);
      G_boid[idx].x  = L_TENTATIVA/3.0f * auxiliar_rand + L_TENTATIVA/2.0f;
      auxiliar_rand = curand_uniform(&G_random_state[idx]);
      G_boid[idx].y  = L_TENTATIVA/3.0f * auxiliar_rand + L_TENTATIVA/2.0f;
      G_boid[idx].label = idx;
      G_boid[idx].v0=V1;
    }
    else
    if(idx<N)
    {
      angle = curand_uniform(&G_random_state[idx]);
      angle = TWO_PI * angle;
      G_boid[idx].vx = V2 * cos(angle);
      G_boid[idx].vy = V2 * sin(angle);
      auxiliar_rand = curand_uniform(&G_random_state[idx]);
      G_boid[idx].x  = L_TENTATIVA/3.0f * auxiliar_rand + L_TENTATIVA/2.0f;
      auxiliar_rand = curand_uniform(&G_random_state[idx]);
      G_boid[idx].y  = L_TENTATIVA/3.0f * auxiliar_rand + L_TENTATIVA/2.0f;
      G_boid[idx].label = idx;
      G_boid[idx].v0=V2;
    }
};



void inicializa(particle *boid, int *random_seed)
{
   int i;   
   float angle; 
   for(i=0;i<N1;i++)
     {
	angle = drand48();
	angle *= TWO_PI;
	boid[i].vx = V1 * cos(angle);
	boid[i].vy = V1 * sin(angle);
	boid[i].x  = L_TENTATIVA/3.0f * drand48() + L_TENTATIVA/2.0f;
	boid[i].y  = L_TENTATIVA/3.0f * drand48() + L_TENTATIVA/2.0f;
	boid[i].label = i;
	boid[i].v0 = V1;
     }
   for(i=N1;i<N;i++)
     {
	angle = drand48();
	angle *= TWO_PI;
	boid[i].vx = V2 * cos(angle);
	boid[i].vy = V2 * sin(angle);
	boid[i].x  = L_TENTATIVA/3.0f * drand48() + L_TENTATIVA/2.0f;
	boid[i].y  = L_TENTATIVA/3.0f * drand48() + L_TENTATIVA/2.0f;
	boid[i].label = i;
	boid[i].v0 = V2;
     }
}
