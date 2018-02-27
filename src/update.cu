#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include <unistd.h>
#include <cuda.h>
//#include <cufft.h>

#include <curand.h>
#include <curand_kernel.h>

#include "param.h"
#include "particle.h"

#define TWO_PI 6.283185307f


__global__ void G_update_position(particle *G_boid, float
*G_Fx, float *G_Fy ,float *G_sum_Vx, float *G_sum_Vy, curandState_t
*state, unsigned int * G_random_number)
{
   //Bem melhor sem memoria compartilhada
   int i=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
   //   float x,y;   
   float module_displace;
   float component_x,component_y;
   float angle;
   if (i<N) 
     {
        G_random_number[i] = curand(state) % RAND_MAX;
	angle = TWO_PI * (float) G_random_number[i] /RAND_MAX;
	component_x =  ( G_sum_Vx[i] + G_Fx[i] + ETA * cos(angle));
	component_y =  ( G_sum_Vy[i] + G_Fy[i] + ETA * sin(angle));
	module_displace = sqrt((component_x*component_x)+(component_y*component_y));
	G_boid[i].x = G_boid[i].x + G_boid[i].v0 * component_x/module_displace;
	G_boid[i].y = G_boid[i].y + G_boid[i].v0 * component_y/module_displace;
	//boundary conditions
	if(G_boid[i].x>=L_TENTATIVA)G_boid[i].x=G_boid[i].x-L_TENTATIVA;
	if(G_boid[i].y>=L_TENTATIVA)G_boid[i].y=G_boid[i].y-L_TENTATIVA;
	if(G_boid[i].x<0)G_boid[i].x=G_boid[i].x+L_TENTATIVA;
	if(G_boid[i].y<0)G_boid[i].y=G_boid[i].y+L_TENTATIVA;
	G_boid[i].vx = G_boid[i].v0 * component_x/module_displace;
	G_boid[i].vy = G_boid[i].v0 * component_y/module_displace; 
     }
};


void  update_position(particle *boid, float *fx, float *fy, float *sum_vx, float *sum_vy)
{
   float angle;
   int i;
   float module_displace;
   float componente_x,componente_y;
   for(i=0;i<N;i++)
     {
	angle = drand48() * TWO_PI;
	componente_x =  ( sum_vx[i] + fx[i] + ETA * cos(angle));
	componente_y =  ( sum_vy[i] + fy[i] + ETA * sin(angle));
	module_displace = sqrt((componente_x*componente_x)+(componente_y*componente_y));
	boid[i].x = boid[i].x + boid[i].v0 * componente_x/module_displace;
	boid[i].y = boid[i].y + boid[i].v0 * componente_y/module_displace;
	// condicoes de contorno
	if(boid[i].x>=L_TENTATIVA)boid[i].x=boid[i].x-L_TENTATIVA;
	if(boid[i].y>=L_TENTATIVA)boid[i].y=boid[i].y-L_TENTATIVA;
	if(boid[i].x<0)boid[i].x=boid[i].x+L_TENTATIVA;
	if(boid[i].y<0)boid[i].y=boid[i].y+L_TENTATIVA;
	boid[i].vx = boid[i].v0 * componente_x/module_displace;
	boid[i].vy = boid[i].v0 * componente_y/module_displace; 
     }   
}

