void  save_position(FILE *arquivo, particle *boid, int time_step)
{
   int i;
   //   printf("entrou");
      fprintf(arquivo,"%d \n",time_step);
         for(i=0;i<N;i++)fprintf(arquivo,"%f %f %f %f\n",boid[i].x,boid[i].y,boid[i].vx,boid[i].vy);
};
