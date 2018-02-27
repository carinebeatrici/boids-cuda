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
   
}
;
