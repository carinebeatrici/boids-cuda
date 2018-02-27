/* Modelo de Vicsek 
 * Inicio: 28/07/2016
 * Versão CUDA - 08/09/2016
 * Sem caixas na memoria global
 */ 

#define TWO_PI 6.2830f

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <cuda.h>
#include <cufft.h>
#include <unistd.h>
#include <curand.h>
#include <curand_kernel.h>

#include "param.h"
#include "particle.h"
#include "cuda-param.h"

#include "distance.cuh"
#include "update.cuh"
#include "save_position.cuh"
 

// Cabeçalhos das funções em C (CPU)
float  *Alocate_vector_real (int vector_size);
float  *Free_vector_real (int vector_size, float *vector);
void   save_position(FILE *arquivo, particle *boid, int time_step);

//void  inicializa(particle *boid, int *random_seed);
//void  calculate_distance(particle *boid, float *fx, float *fy, float *sum_vx, float *sum_vy);
//void  update_position(particle *boid, float *fx, float *fy, float *sum_vx, float *sum_vy);

//Cabeçalh das funçõ em CUDA (GPU)
__global__ void G_initialize(particle *G_boid,  curandState_t *G_random_state, unsigned int *G_random_number);
__global__ void G_initialize_random_generator(unsigned int seed, curandState_t *state);
//__global__ void G_update_position(particle *G_boid, float *G_Fx, float *G_Fy ,float *G_sum_Vx, float *G_sum_Vy, curandState_t *state, unsigned int *G_random_number);
//__global__ void G_calculate_distance(particle *G_boid, float *G_Fx, float *G_Fy ,float *G_sum_Vx, float *G_sum_Vy);



int main (void)
{
//   printf("Começo do programa  \n");
//   float *fx,*fy,*sum_vx,*sum_vy; // coordinates, velocities and forces projections
   float exec_time_gpu=0.0,exec_time_gpu_total=0.0,exec_time_cpu=0.0;
   float exec_time_total=0.0;
   clock_t time_cpu_init,time_cpu_end;
   cudaEvent_t start,stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   time_cpu_init=clock();
   float *G_fx,*G_fy,*G_sum_vx,*G_sum_vy; // coordinates, velocities and forces projections on GPU
   particle *boid;   //particulas em CPU para salvar os arquivos
   particle *G_boid; //particulas em GPU para calculos
   curandState_t *G_random_state; //sementes aleatorias
   unsigned int *G_random_number;  // vetores com numeros aleatorios em GPU
   long int time_step=0;
   FILE *arquivo_saida;
   arquivo_saida=fopen("data/posicoes.dat","w");
   // Allocating stuff   
//   printf(" Alocando coisas...................................");
   boid = (particle*) malloc (N * sizeof(particle)); // Alocando vetor de estruturas em CPU
   cudaMalloc((void**) &G_boid, (N * sizeof(particle))); //N random numbers each time step
   cudaMalloc((void**) &G_random_state, N * sizeof(curandState_t)); //N random numbers each time step
   cudaMalloc((void**) &G_random_number, N * sizeof(unsigned int)); //N random numbers each time step   
   cudaMalloc((void**) &G_fx, N * sizeof(float)); 
   cudaMalloc((void**) &G_fy, N * sizeof(float)); 
   cudaMalloc((void**) &G_sum_vx, N * sizeof(float));
   cudaMalloc((void**) &G_sum_vy, N * sizeof(float));
//   printf("OK  \n");
//   fx     = Alocate_vector_real (N);
//   fy     = Alocate_vector_real (N);
//   sum_vx = Alocate_vector_real (N);
//   sum_vy = Alocate_vector_real (N);

//   inicializa(boid);
//   printf(" Antes inicializar random generator................");
   time_cpu_end = clock();
   exec_time_cpu = time_cpu_end-time_cpu_init;
   time_cpu_init=clock();
   cudaEventRecord(start,0);
   G_initialize_random_generator<<<BLOCKS,THREADS>>>(time(0),G_random_state);
//   printf("OK  \n");
// Inicia sementes aleatorias em GPU
//   printf("Inicializando boids................................");
   G_initialize<<<BLOCKS,THREADS>>>(G_boid,G_random_state,G_random_number);
//   printf("OK  \n");
//   printf("Copiando dados.....................................");
   cudaMemcpy(boid,G_boid,(N*sizeof(particle)),cudaMemcpyDeviceToHost);
   cudaEventRecord(stop,0);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&exec_time_gpu,start,stop);
   exec_time_gpu_total+=exec_time_gpu;
   time_cpu_init=clock();
//   printf("OK  \n");
//   printf("x: %f y: %f vx: %f vy: %f  \n",boid[0].x,boid[0].y,boid[0].vx,boid[0].vy);
//   printf("x: %f y: %f vx: %f vy: %f  \n",boid[N-1].x,boid[N-1].y,boid[N-1].vx,boid[N-1].vy);
//   printf("Salvando Particulas................................");
   save_position(arquivo_saida, boid, time_step);
//   printf("OK  \n");
//   printf("Stuff alocated \n");
   // Time loop
   while (time_step < TIME_FINAL)
     {
//	printf("Loop time...........................................,%d  \n",time_step);
	time_step++;
//	calculate_distance(boid,fx,fy,sum_vx,sum_vy);
	time_cpu_end=clock();
	exec_time_cpu += time_cpu_end-time_cpu_init;
	cudaEventRecord(start,0);
        G_calculate_distance<<<BLOCKS,THREADS>>>(G_boid,G_fx,G_fy,G_sum_vx,G_sum_vy);
//	update_position(boid,fx,fy,sum_vx,sum_vy);
        G_update_position<<<BLOCKS,THREADS>>>(G_boid,G_fx,G_fy,G_sum_vx,G_sum_vy,G_random_state,G_random_number);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&exec_time_gpu,start,stop); 
	exec_time_gpu_total+=exec_time_gpu;
	time_cpu_init=clock();
	if(time_step%1000==0)
	  {
	    cudaMemcpy(boid,G_boid,N*sizeof(particle),cudaMemcpyDeviceToHost);
	    save_position(arquivo_saida, boid, time_step);
	  }
     }
   // LIBERATING MEMORY
//   printf("Stuff dealocated \n");
//   fx = Free_vector_real (N, fx);
//   fy = Free_vector_real (N, fy);
//   sum_vx = Free_vector_real (N, sum_vx);
//   sum_vy = Free_vector_real (N, sum_vy);
   cudaFree(G_boid);
   cudaFree(G_random_state);
   cudaFree(G_random_number);
   cudaFree(G_fx);
   cudaFree(G_fy);
   cudaFree(G_sum_vx);
   cudaFree(G_sum_vy);
   fclose(arquivo_saida);
   time_cpu_end=clock();
   exec_time_cpu += time_cpu_end-time_cpu_init;
   exec_time_cpu /= (float)CLOCKS_PER_SEC;
   exec_time_gpu_total /= 1000.0;
   exec_time_total= exec_time_cpu+ exec_time_gpu_total;
   printf("   \n");
   printf("tempo total: %f s numero de particulas: %d \n tempo de CPU: %f s tempo de GPU: %f s\n",
	  exec_time_total, N, exec_time_cpu, exec_time_gpu_total );
   printf("Threads: %d Blocks: %d \n",THREADS,BLOCKS);
//   printf("Num of loads: %d \n",NUM_LOADS);
   return 0;
};




/*
void  save_position(FILE *arquivo, particle *boid, int time_step)
{
   int i;
//   printf("entrou");
   fprintf(arquivo,"%d \n",time_step);
   for(i=0;i<N;i++)fprintf(arquivo,"%f %f %f %f\n",boid[i].x,boid[i].y,boid[i].vx,boid[i].vy);
};
*/

/*
__global__ void G_calculate_distance(particle *G_boid, float *G_Fx, float *G_Fy ,float *G_sum_Vx, float *G_sum_Vy)
{
   //---------------locais--------------------------------------------------
   int i,j;
   float dx,dy,L_TENTATIVAh=L_TENTATIVA*0.5f,distance;
   float r_max2=R_MAX*R_MAX;
   float auxiliar;
   //O loop inplicito vai ser na variavel i
   //Cada thread vai calcular as forcas para uma particula
   i = (blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
   //------------------todo mundo com todo  mundo----------------------------
   //------------------dentro da memoria global-----------------------------
   if(i<N)
     {
	G_Fx[i]=0.0f;
	G_Fy[i]=0.0f;
	G_sum_Vx[i]=0.0f;
	G_sum_Vy[i]=0.0f;
	for(j=0;j<N;j++)
	  {
	     if(i!=j)
	       {
		  dx=G_boid[i].x-G_boid[j].x;
		  dy=G_boid[i].y-G_boid[j].y;
		  if(dx>L_TENTATIVAh)dx=dx-L_TENTATIVA;else if(dx<-L_TENTATIVAh) dx=L_TENTATIVA+dx;
		  if(dy>L_TENTATIVAh)dy=dy-L_TENTATIVA;else if(dy<-L_TENTATIVAh) dy=L_TENTATIVA+dy;
		  distance=dx*dx+dy*dy;
		  //check if particles are interacting
		  if(distance<r_max2)
		    {
		       //align
		       if(j<N1&&i<N1)
			 {
			    G_sum_Vx[i]+= ALPHA11 * G_boid[j].vx;
			    G_sum_Vy[i]+= ALPHA11 * G_boid[j].vy;
			 }
		       else
			 if(j>N1&&i>N1)
			   {
			      G_sum_Vx[i]+= ALPHA22 * G_boid[j].vx;
			      G_sum_Vy[i]+= ALPHA22 * G_boid[j].vy;
			   }
		       else
			 {
			    G_sum_Vx[i]+= ALPHA12 * G_boid[j].vx;
			    G_sum_Vy[i]+= ALPHA12 * G_boid[j].vy;
			 }
		       distance=sqrt(distance);
		       //hard core replusion
		       if (distance<=R_CORE)
			 { //hard core replusion
			    auxiliar=FORCA_CORE/distance;
			    G_Fx[i]+=dx*auxiliar;
			    G_Fy[i]+=dy*auxiliar;
			 }
		       else
			 {
			    //valid force range
			    if(j<N1&&i<N1)
			      {
				 G_Fx[i]+=BETA11*dx*(1.0f-distance/R_EQ);
				 G_Fy[i]+=BETA11*dy*(1.0f-distance/R_EQ);
			      }
			    else
			      {
				 if(i>=N1&&j>=N1)
				   {
				      G_Fx[i]+=BETA22*dx*(1.0f-distance/R_EQ);
				      G_Fy[i]+=BETA22*dy*(1.0f-distance/R_EQ);
				   }
				 else
				   {
				      if((i>=N1&&j<N1)||(i<N1&&j>=N1))
					{
					   G_Fx[i]+=BETA12*dx*(1.0f-distance/R_EQ);
					   G_Fy[i]+=BETA12*dy*(1.0f-distance/R_EQ);
					}
				   }
			      }
			 }
		    }
	       }
	  }
     }
};
*/   

/*
__global__ void G_update_position(particle *G_boid, float *G_Fx, float *G_Fy ,float *G_sum_Vx, 
				  float *G_sum_Vy, curandState_t *state, unsigned int * G_random_number)
{
   int i=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
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
  */ 