// Calculo das distancias e forças entre as particulas dentro da GPU

#include <math.h>
//#include <unistd.h>
#include <cuda.h>
//#include <cufft.h>
#include <math.h>

#include "param.h"
#include "particle.h"

	 

__global__ void G_distance(particle *G_boid, float *G_Fx, float *G_Fy ,float *G_sum_Vx, float *G_sum_Vy)
{
   //---------------locais--------------------------------------------------
   int i,j;
   float dx,dy,L_TENTATIVAh=L_TENTATIVA*0.5f,distance;
   float r_max2=R_MAX*R_MAX;
   float auxiliar;
   //O loop inplicito vai ser na variavel i
   //Cada thread vai calcular as forcas para uma particula
   i = (blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
   //------------------todo mundo com todo mundo----------------------------
   //------------------dentro da memoria global-----------------------------
   G_Fx[i]=0.0f;
   G_Fy[i]=0.0f;
   G_sum_Vx[i]=0.0f;
   G_sum_Vy[i]=0.0f;
   for(j=0;j<N;j++)
     {     
	if(i!=j)
	  {
	     dx=G_boid[j].x-G_boid[i].x;
	     dy=G_boid[j].y-G_boid[i].y;
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
		       auxiliar=-FORCA_CORE/distance;
		       G_Fx[i]+=dx*auxiliar;
		       G_Fy[i]+=dy*auxiliar;
		    }
		  else
		    {
		       //valid force range
		       if(j<N1&&i<N1)
			 {
			    G_Fx[i]-=BETA11*dx*(1.0f/distance-FORCA_CORE);
			    G_Fy[i]-=BETA11*dy*(1.0f/distance-FORCA_CORE);
			 }
		       else
			 {
			    if(i>=N1&&j>=N1)
			      {
				 G_Fx[i]-=BETA22*dx*(1.0f/distance-FORCA_CORE);
				 G_Fy[i]-=BETA22*dy*(1.0f/distance-FORCA_CORE);
			      }
			    else
			      {
				 if((i>=N1&&j<N1)||(i<N1&&j>=N1))
				   {
				      G_Fx[i]-=BETA12*dx*(1.0f/distance-FORCA_CORE);
				      G_Fy[i]-=BETA12*dy*(1.0f/distance-FORCA_CORE);
				   }	 
			      }  
			 }
		    }
	       }
	  }
     } 
};


void calculate_distance(particle *boid , float *fx, float *fy, float *sum_vx, float *sum_vy)
{
   int i,j;
   float forca;
   float dx,dy;
   float distance;
   // type 1 and type 1
   for(i=0;i<N1;i++)
     {
	fx[i]=0.00f;
	fy[i]=0.00f;
	sum_vx[i]=0.00f;
	sum_vy[i]=0.00f;
	for(j=0;j<N1;j++)
	  {
	     dx=boid[i].x-boid[j].x;
	     dy=boid[i].y-boid[j].y;
	     if(dx>L_TENTATIVA/2.0f)dx=dx-L_TENTATIVA;
	     if(dy>L_TENTATIVA/2.0f)dy=dy-L_TENTATIVA;
	     if(dx<-L_TENTATIVA/2.0f)dx=dx+L_TENTATIVA;
	     if(dy<-L_TENTATIVA/2.0f)dy=dy+L_TENTATIVA;
	     distance=sqrt((dx*dx)+(dy*dy));
	     if(distance<R_CORE)
	       {
		  forca = FORCA_CORE; 
		  if(distance>0.00010f)
		    {
		       sum_vx[i]+= ALPHA11 * boid[j].vx;
		       sum_vy[i]+= ALPHA11 * boid[j].vy;
		    }		  
	       }
	     else 
	       {
		  if(distance<R_MAX)
		    {
		       forca = BETA11 * (1.00f - distance/R_EQ);
		       sum_vx[i]+= ALPHA11 * boid[j].vx;
		       sum_vy[i]+= ALPHA11 * boid[j].vy;
		    }
		  else
		    {
		       forca = 0.00f;
		    }
	       }
	     if(distance>0.00010f)
	       {
		  fx[i] += (forca * dx) / distance;
		  fy[i] += (forca * dy) / distance;
	       }
	  }
     }
   // type 1 and type 2
   for(i=0;i<N1;i++)
     {
	for(j=N1;j<N;j++)
	  {
	     dx=boid[i].x-boid[j].x;
	     dy=boid[i].y-boid[j].y;
	     if(dx>L_TENTATIVA/2.0f)dx=dx-L_TENTATIVA;
	     if(dy>L_TENTATIVA/2.0f)dy=dy-L_TENTATIVA;
	     if(dx<-L_TENTATIVA/2.0f)dx=dx+L_TENTATIVA;
	     if(dy<-L_TENTATIVA/2.0f)dy=dy+L_TENTATIVA;
	     distance=sqrt((dx*dx)+(dy*dy));
	     if(distance<R_CORE)
	       {
		  forca = FORCA_CORE; 
		  if(distance>0.00010f)
		    {
		       sum_vx[i]+= ALPHA12 * boid[j].vx;
		       sum_vy[i]+= ALPHA12 * boid[j].vy;
		    }		  
	       }
	     else 
	       {
		  if(distance<R_MAX)
		    {
		       forca = BETA12 * (1.00f - distance/R_EQ);
		       sum_vx[i]+= ALPHA12 * boid[j].vx;
		       sum_vy[i]+= ALPHA12 * boid[j].vy;
		    }
		  else
		    {
		       forca = 0.00f;
		    }
	       }
	     if(distance>0.00010f)
	       {
		  fx[i] += (forca * dx) / distance;
		  fy[i] += (forca * dy) / distance;
	       }
	  }
     }
   // type 2 and type 1
   for(i=N1;i<N;i++)
     {
	fx[i]=0.00f;
	fy[i]=0.00f;
	sum_vx[i]=0.00f;
	sum_vy[i]=0.00f;
	for(j=0;j<N1;j++)
	  {
	     dx=boid[i].x-boid[j].x;
	     dy=boid[i].y-boid[j].y;
	     if(dx>L_TENTATIVA/2.0f)dx=dx-L_TENTATIVA;
	     if(dy>L_TENTATIVA/2.0f)dy=dy-L_TENTATIVA;
	     if(dx<-L_TENTATIVA/2.0f)dx=dx+L_TENTATIVA;
	     if(dy<-L_TENTATIVA/2.0f)dy=dy+L_TENTATIVA;
	     distance=sqrt((dx*dx)+(dy*dy));
	     if(distance<R_CORE)
	       {
		  forca = FORCA_CORE; 
		  if(distance>0.00010f)
		    {
		       sum_vx[i]+= ALPHA12 * boid[j].vx;
		       sum_vy[i]+= ALPHA12 * boid[j].vy;
		    }		  
	       }
	     else 
	       {
		  if(distance<R_MAX)
		    {
		       forca = BETA12 * (1.00f - distance/R_EQ);
		       sum_vx[i]+= ALPHA12 * boid[j].vx;
		       sum_vy[i]+= ALPHA12 * boid[j].vy;
		    }
		  else
		    {
		       forca = 0.00f;
		    }
	       }
	     if(distance>0.00010f)
	       {
		  fx[i] += (forca * dx) / distance;
		  fy[i] += (forca * dy) / distance;
	       }
	  }
     }
   // type 2 and type 2 
   for(i=N1;i<N;i++)
     {
	for(j=N1;j<N;j++)
	  {
	     dx=boid[i].x-boid[j].x;
	     dy=boid[i].y-boid[j].y;
	     if(dx>L_TENTATIVA/2.0f)dx=dx-L_TENTATIVA;
	     if(dy>L_TENTATIVA/2.0f)dy=dy-L_TENTATIVA;
	     if(dx<-L_TENTATIVA/2.0f)dx=dx+L_TENTATIVA;
	     if(dy<-L_TENTATIVA/2.0f)dy=dy+L_TENTATIVA;
	     distance=sqrt((dx*dx)+(dy*dy));
	     if(distance<R_CORE)
	       {
		  forca = FORCA_CORE; 
		  if(distance>0.00010f)
		    {
		       sum_vx[i]+= ALPHA22 * boid[j].vx;
		       sum_vy[i]+= ALPHA22 * boid[j].vy;
		    }		  
	       }
	     else 
	       {
		  if(distance<R_MAX)
		    {
		       forca = BETA22 * (1.00f - distance/R_EQ);
		       sum_vx[i]+= ALPHA22 * boid[j].vx;
		       sum_vy[i]+= ALPHA22 * boid[j].vy;
		    }
		  else
		    {
		       forca = 0.00f;
		    }
	       }
	     if(distance>0.00010f)
	       {
		  fx[i] += (forca * dx) / distance;
		  fy[i] += (forca * dy) / distance;
	       }
	  }
     }
}
