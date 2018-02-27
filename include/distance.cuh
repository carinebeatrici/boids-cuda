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
   //------------------todo mundo com todo mundo----------------------------
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
		  if(dx>L_TENTATIVAh)dx=dx-L_TENTATIVA;else
		    if(dx<-L_TENTATIVAh) dx=L_TENTATIVA+dx;
                  if(dy>L_TENTATIVAh)dy=dy-L_TENTATIVA;else
		    if(dy<-L_TENTATIVAh) dy=L_TENTATIVA+dy;
                  distance=dx*dx+dy*dy;
		  //check if particles are interacting
		  if(distance<r_max2)
		    {
		       //align
		       if(j<N1&&i<N1)			                                {
			  G_sum_Vx[i]+= ALPHA11 * G_boid[j].vx;
			  G_sum_Vy[i]+=ALPHA11 * G_boid[j].vy;
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
			 {
			    //hard core replusion
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

