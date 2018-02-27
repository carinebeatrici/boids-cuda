/*
 *  tentativas para otimizacao utilizando a memoria compartilhada
 *  calcuco do gama feito ok mas nao ficou mais rapido :(
 *  Ficou mais lento precisamos otimizar mesmo é a update
 *  Ficou mais rápido ajustando o número de blocos e threads
 *  ~2.4x mais rapido de 5.91s para 2.5s
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <unistd.h>
#include <cuda.h>
#include <cufft.h>
#include "paramfile.h"
#include "fileutils.h"
#include "stringutils.h"

#include "particles.h"

using namespace std;

#define N_tot 6000
//===========================================================================
//===========================================================================
//===========================================================================
// CUDA  kernel&device functions
//---------------------------------------------------------------------------

__global__ void G_setrconst(int N,REAL *a,REAL c) {
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    if(idx<N) {a[idx]=c;}
};

__global__ void G_seticonst(int N,int *a,int c) {
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    if(idx<N) {a[idx]=c;}
};

//---------------------------------------------------------------------------
__device__ int Gran0(int rs)
{
	int k;
	k=rs/IQ0;
	rs=IA0*(rs-k*IQ0)-IR0*k;
	if (rs < 0) rs += IM0;
	return rs;
};

__global__ void G_addrand(int N,REAL *A,int *rs,REAL T)
{
//calculo dos numeros aleatorios 
        int i=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
	int rn;
	if(i<N)
	{
		rn=rs[i];
		rn=Gran0(rn);
		A[i]+=T*(2.0*AM0*rn-1.0);
		rs[i]=rn;
	}
};
//---------------------------------------------------------------------------

#define FMAX 1000.0
#define b11 0.6
#define b12 0.12 
#define b22 0.1
#define a11 0.0 
#define a12 0.0 
#define a22 0.0 
#define slope11 2.5
#define slope22 2.5
#define slope12 (slope11+slope22)/2.0
#define constforca11 slope11*0.4
#define constforca12 slope12*0.4
#define constforca22 slope22*0.4

__global__ void G_update(int N,part *Gp,REAL alpha,REAL beta, REAL eta,REAL a2,REAL f0,REAL ra2,REAL L,int *rs,int Nslow)
{
    int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
    int i,rn1,rn2;
    __shared__ REAL Sx[N_tot];
    __shared__ REAL Sy[N_tot];
    if(threadIdx.x==0)
     {	
	for(i=0;i<N_tot;i++)
	  {	
	     Sx[i]=Gp[i].x;
	     Sy[i]=Gp[i].y;
	  }
     }
    __syncthreads();
    REAL dx,dy,Lh,r,r2,Fx,Fy,b,Vx,Vy,Tx,Ty;
    if(idx<N_tot){
      Lh=L*0.5;
      Fx=0.0;
      Fy=0.0;
      Vx=0.0;
      Vy=0.0;
      for (i=0; i<N; i++) {
        if(i!=idx)
	{   // calculate distance taking PBC into account
	  dx=Sx[idx]-Sx[i];
	  dy=Sy[idx]-Sy[i];
	  // Periodic boundary conditions
	  if(dx>Lh)dx=dx-L;   else if(dx<-Lh)dx=L+dx;
	  if(dy>Lh)dy=dy-L;   else if(dy<-Lh)dy=L+dy;
	  r2=dx*dx+dy*dy;
	  //check if particles are interacting
	  if(r2<ra2)
	  //if(r2<0.3025)
	  {  
	    r=sqrt(r2);
	    if (r2<=a2) 
	    { //hard core replusion
	      b=FMAX/r;
	      Fx+=dx*b;
	      Fy+=dy*b;
	    }
	    else
	    { //valid force range
	       if(idx<Nslow&&i<Nslow)
	       {
	       Fx+=b11*dx*((constforca11/r)-slope11);
	       Fy+=b11*dy*((constforca11/r)-slope11);
	       }
	       else
	       {
	         if(idx>=Nslow&&i>=Nslow)
	         {
		    Fx+=b22*dx*((constforca22/r)-slope22);
		    Fy+=b22*dy*((constforca22/r)-slope22);
                 }
		 else
		 {
		   if((idx>=Nslow&&i<Nslow)||(idx<Nslow&&i>=Nslow))
	           {
		     Fx+=b12*dx*((constforca12/r)-slope12);
		     Fy+=b12*dy*((constforca12/r)-slope12);
		   }
		 }		                                     
	      }
	   }      
	   //sum of neighbor velocities
	   if(idx<Nslow&&i<Nslow)
	   {
	      Vx+=Gp[i].vx*a11;
              Vy+=Gp[i].vy*a11;
	   }
	   if(idx>=Nslow&&i>=Nslow)
	   {
	      Vx+=Gp[i].vx*a22;
	      Vy+=Gp[i].vy*a22;
	   }
	   if((idx>=Nslow&&i<Nslow)||(idx<Nslow&&i>=Nslow))
	   {
	      Vx+=Gp[i].vx*a12;
	      Vy+=Gp[i].vy*a12;
	   }   
	   }
	}

   }
      //random angle
      rn1=Gran0(rs[idx]);
      rs[idx]=rn1;
      b=TWOPI*(AM0*rn1)-PI;
      rn2=Gran0(rs[idx]);
      rs[idx]=rn2;
      Tx=Vx+Fx+eta*cos(b);
      Ty=Vy+Fy+eta*sin(b);
      Gp[idx].theta=atan2(Ty,Tx);      
    } 
};


__global__ void G_timestep(int N,part *Gp,REAL dt,REAL L) {
//Bem melhor sem memoria compartilhada
   int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
   REAL vx,vy,theta,v0,x,y;   
   if (idx<N) {
      //update velocity
      theta=Gp[idx].theta;
      v0=Gp[idx].v0;
      vx=v0*cos(theta);
      vy=v0*sin(theta);
      Gp[idx].vx=vx;
      Gp[idx].vy=vy;
      //update coordinate
      x=Gp[idx].x;
      y=Gp[idx].y;
      x+=dt*vx;while(x<0) x+=L;while(x>=L) x-=L;
      y+=dt*vy;while(y<0) y+=L;while(y>=L) y-=L;
      Gp[idx].x=x;
      Gp[idx].y=y;
      }
};

__global__ void G_gama(int N,part *Gp, REAL ra2,REAL Nslow)
{
__shared__ int Sna[N_tot];
__shared__ int Snv[N_tot];   
int idx=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
REAL dx,dy,r2;
int i;
if(idx<N){
   Sna[idx]=0;
   Snv[idx]=0;
}
if(idx<N){
   for(i=0;i<N;i++){
      if(i!=idx){
         dx=Gp[idx].x-Gp[i].x;
	 dy=Gp[idx].y-Gp[i].y;
	 r2=dx*dx+dy*dy;
	 //if(r2<ra2){
	 if(r2<0.3025){
	    if(i<Nslow)Snv[idx]++;
	    else Sna[idx]++;
	 }
      }
   }
}
//__syncthreads();
if(idx<N)
     {	
	Gp[idx].na=Sna[idx];
	Gp[idx].nv=Snv[idx];
     }
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//-------parametros----------------------------------------------------------
particles::particles()
{  
    idum=2354;
    iyy=0;
    ir=new int[100];  iff=0;
    simtype=0;
    N=N_tot;   //total number
    Nslow=3000; //number of slow particles
    Nt=100000000; //number of timesteps
    Nto=1000; //output interval
    dt=1.0;
    vslow=0.007; //slow velocity
    vfast=0.007; //fast velocity
    alpha=0.0;
    beta=0.55;
    eta=1.0;
//    L=(2.*ceil(pow(N,0.5)*0.17))+4.; //system size
    L=19;
    a=0.2;  //particle size
    ra=0.55; //force radius
    f0=2.5; //force slope
    rs=123456; //random seed
    testecontinua=0;
    Gdev=0;
    outfn="part";
    Gp=NULL;
    GRs=NULL;
    p=NULL;
};


particles::~particles()
{
    if(GRs!=NULL) cudaFree(GRs);
    if(Gp!=NULL) cudaFree(Gp);
    if(p!=NULL) delete[] p;
    delete[] ir;
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double particles::ran2(int iseed)
{
	//       based on ran2 from Numerical Recipes page 197 (fortran)
	//       meant to emulate the unix function rand
	//       positive iseed reinitializes, zero iseed - same string
	// ir - class-array int[100], ma,ia,ic,rm - constants
	// idum,iff,iyy - class-var (iff=0 before first run)
	double rand;
	int j;
	if((iseed>0) || (iff==0))
	{iff=1;
		idum=-iseed-3;
		idum=(ic-idum)%ma;       //-> idum < ma
		for(j=0;j<97;j++)  {
			idum=(ia*idum+ic)%ma;      //-> idum < ma
			ir[j]=idum;
		}
		idum=(ia*idum+ic)%ma;         //-> idum < ma
		iyy=idum;                     //-> iyy <ma
	}
	j=(97*iyy)/ma;                  //-> i < 97
	iyy=ir[j];
	rand=iyy*rm;
	idum=(ia*idum+ic)%ma;
	ir[j]=idum;
	if(rand==0.0) rand=1e-14; //double !!
	return rand;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void particles::simulate()
{
    int frame,i;
    double t;
    long int tnf=1,n,tempocontinuacao;
    REAL gama=0;
    FILE *arq1;
    FILE *arq2;
    FILE *arq3;
    frame=0;
    t=0.0;
//    abrindo o arquivo com informacoes do estado final do sistema interrompido
//  tudo ainda em cpu... 
   if(testecontinua==1)
     {
	arq3=fopen("estado-final","r");
	fscanf(arq3,"%d %d %li %d %d %f %f %f %f %f %f",&N,&Nslow,&n,&tnf,&frame,&alpha,&beta,&eta,&a,&ra,&L);	
	printf("%d %d %d %li %d %f %f %f %f %f %f\n",N,Nslow,n,tnf,frame,alpha,beta,eta,a,ra,L);
	//alocando vetor de p na CPU
	p=new part[N];
	//alocando Gp na GPU
	memP=N*sizeof(part);
	cudaMalloc((void**)&Gp,memP);
	printf("%f %f %f\n",f0,vslow,vfast);
	printf("N= %d, Nslow=%d\n",N,Nslow);
	for(i=0;i<N;i++)
	  {
	     fscanf(arq3,"%f %f %f %f",&p[i].x,&p[i].y,&p[i].vx,&p[i].vy);
	  }
	for(i=0;i<Nslow;i++)p[i].v0=vslow;
	for(i=Nslow;i<N;i++)p[i].v0=vfast;
	t=n*dt;
	arq2=fopen("posicoes.dat","a");
	// copia os dados de p* para a GPU
	cudaMemcpy(Gp,p,memP,cudaMemcpyHostToDevice);//(destino,fonte,memoria,cudaMemcpyHostToDevice ou cudaMemcpyDeviceToHost)	
	tempocontinuacao=n;
     }
   else
     {	
	arq1=fopen("dados.dat","w");
	fclose(arq1);
	arq2=fopen("posicoes.dat","w");
	fclose(arq2);
	output(frame,gama,t); //intial state, data is still on CPU
	tempocontinuacao=1;
     }
//    printf("alfa %f, beta %f\n",alpha,beta);
    for(n=tempocontinuacao;n<=Nt;n++)
    {        
        G_update<<<GRID,BLOCK>>>(N,Gp,alpha,beta,eta,a,f0,ra,L,GRs,Nslow); //update particle velocities
        //cudaMemcpy(p,Gp,memP,cudaMemcpyDeviceToHost); //copia da GPU para CPU 
        G_timestep<<<GRID,BLOCK>>>(N,Gp,dt,L); //integrate equation of motion
        t+=dt;
	if((n==1)||(n%200)==0)
        {	   
	   arq1=fopen("dados.dat","a");    
	   G_gama<<<GRID,BLOCK>>>(N,Gp,ra,Nslow); //calculo dos numeros de vizinhos na e nv
	   cudaMemcpy(p,Gp,memP,cudaMemcpyDeviceToHost); //copia da GPU para CPU 
	   gama=0;
	   for(i=0;i<Nslow;i++){if((p[i].na+p[i].nv)!=0)gama=gama+(p[i].na/(p[i].na+p[i].nv));}
	   gama/=Nslow;
	   outputgama(gama,n,arq1);
	   printf("n=%d, gama=%f\n",n,gama);
	   fflush(stdout);
	   fclose(arq1);
	   if((n%Nto==0)||(n==1))
	     {
		arq3=fopen("estado-final","w");
		outputfinal(N,Nslow,n,tnf,frame,alpha,beta,eta,a,ra,L,arq3);
		fclose(arq3);
	     }
        }	
	if(n>=tnf){
	   tnf=tnf+pow(tnf,0.8);
	   cudaMemcpy(p,Gp,memP,cudaMemcpyDeviceToHost);
	   arq2=fopen("posicoes.dat","a");
	   outputposition(n,arq2);
	   fclose(arq2);
	   fflush(stdout);
	}							
        if((n%Nto)==0)
        {
	    frame++;
            cudaMemcpy(p,Gp,memP,cudaMemcpyDeviceToHost);
	    output(frame,gama,t);
        }
    }
//fclose(arq1);
};

//---------------------------------------------------------------------------
int particles::output(int fno, REAL gama, REAL t)
{
//    int i;
//    fHandle f;
//    string s;

//    f=FileCreate(outfn+IntToStr(fno)+".dat");
//    s="#t= "+FloatToStr(t)+", N= "+IntToStr(N)+", L= "+FloatToStr(L)+", gama= "+FloatToStr(gama)+" a= "+FloatToStr(a)+"\n";
//    FileWrite(f,s.c_str(),s.length());
    
//    for(i=0;i<N;i++)
//    {
//        s=IntToStr(i+1)+"\t"+FloatToStr(p[i].x)+"\t"+FloatToStr(p[i].y)+"\t"+FloatToStr(p[i].v0)+"\t"+FloatToStr(p[i].theta)+"\n";
//    	FileWrite(f,s.c_str(),s.length());
//    }
    
//    FileClose(f);
    return 0;
};

//--------------------------------------------------------------------------
int particles::outputgama(REAL gama, int passo, FILE *arq1)
{
    fprintf(arq1,"%d %f\n",passo,gama); fflush(stdout);
    return 0;
};

//--------------------------------------------------------------------------
int particles::outputfinal(int N,int Nslow,long int n,long int tnf, int frame, float alpha, float beta, float eta, float a, float ra, float L, FILE *arq3)
{
   int i;
   // imprimindo parametros do sistema
   // write(30,*)l,tnf,q,a,raz
   // write(30,"(4(f10.6,1x))")x,y
   //printf("escreveu estado-final\n");
   fprintf(arq3,"%d %d %li %li %d %f %f %f %f %f %f \n",N,Nslow ,n,tnf,frame,alpha, beta,eta,a,ra,L); 
   //escreve posicoes e angulo
   for(i=0;i<N;i++)fprintf(arq3,"%f %f %f %f\n",p[i].x,p[i].y,p[i].vx,p[i].vy);
   fflush(stdout);
   return 0;
};


//--------------------------------------------------------------------------
int particles::outputposition(int n,FILE *arq2)
{
    int i;
    fprintf(arq2,"%d\n",n);
    for(i=0;i<Nslow;i++)fprintf(arq2,"%f %f %f %f\n",p[i].x,p[i].y,p[i].vx,p[i].vy);
    fprintf(arq2,"inicio ciano\n");
    for(i=Nslow;i<N;i++)fprintf(arq2,"%f %f %f %f\n",p[i].x,p[i].y,p[i].vx,p[i].vy);
    return 0;
}
			
//---------------------------------------------------------------------------
void particles::initparams(int n, char* argv[])
{//overwrite deault built-in parameters set in constructor
	int m; 
	paramfilereader *pf=new paramfilereader();
	m=0;
	if(n>1) m=pf->cmdlineopenread(n,argv);
	if(m<1) m=pf->openread("part.ini");
	if(m>0)
	{  
	   simtype=pf->getint("simtype",simtype);
	   Gdev=pf->getint("Gdev",Gdev);
	   testecontinua=pf->getint("c",testecontinua);
	   N=pf->getint("N",N);
	   Nslow=pf->getint("Nslow",Nslow);
	   Nt=pf->getint("Nt",Nt);
	   Nto=pf->getint("Nto",Nto);
	   dt=pf->getdouble("dt",dt);
	   L=pf->getdouble("L",L);
	   a=pf->getdouble("a",a);
	   alpha=pf->getdouble("alpha",alpha);
	   beta=pf->getdouble("beta",beta);
	   eta=pf->getdouble("eta",eta);
	   outfn=pf->getstring("output");
	   if(outfn=="") outfn="part";		

	}
	delete pf;
	
};

//---------------------------------------------------------------------------

int particles::initialize()
{
    double x;
    REAL angle;//,angulo, raio;
//    REAL raiosorteio;
    int i,j,k;
    cudaError_t err;
    size_t fmem,tmem;
    int *iarr;

    cudaThreadExit();
    cudaSetDevice(Gdev);
    cuMemGetInfo(&fmem,&tmem);
    printf("GPU memory before allocation free: %u, total: %u\n",fmem,tmem);
		
    BLOCK=dim3(MAXT,1);
    // blocs 512, 1024 , tentando 8000
    i=8; //x blocks, limited to 256^2 !!! because of the GPU and the cuda version
    //testando: i=8000 antes de mexer hoje (13/03/2014)
    k=i*BLOCK.x;
    j=(N+k-1)/k; //y blocks
    GRID=dim3(i,j);
    
    i=MAXT;
    printf("GPU execution layout:\n - threads: %d\n - system block: %dx%d\n - state block %d\n",i,GRID.x,GRID.y,BLOCK.x);
	
    memP=N*sizeof(part);
    memRs=N*sizeof(int);

    cudaMalloc((void**)&Gp,memP);
    cudaMalloc((void**)&GRs,memRs);
	
    ran2(rs);	
    iarr=new int[N];
	 
    for(i=0;i<N;i++) 
    { 
        x=ran2()*0xFFFFFFF+127983.0;
    	j=((int) x)^MASK0;
    	iarr[i]= j;
    }
    cudaMemcpy(GRs, iarr,memRs,cudaMemcpyHostToDevice);
    delete[] iarr;
	
    err=cudaGetLastError();
    if(err!=cudaSuccess)                                 
    printf("CUDA error [%d] (alloc) : %s\n",err,cudaGetErrorString(err));
	
    cuMemGetInfo(&fmem,&tmem);
    printf("GPU memory after allocation free: %u, total: %u\n",fmem,tmem);
    p=new part[N];	
    if(testecontinua==0)
      {   	
	//sorteio inicial    
	for(i=0;i<N;i++)
	  {
	     if(i<Nslow) p[i].v0=vslow;
	     else p[i].v0=vfast;
	   //  raiosorteio=(L-5.)/2.;
	     // sorteio inicial redondo
	    // angulo=ran2()*TWOPI;
	    // raio=pow(ran2(),0.5);
	    // p[i].x=L/2.+raiosorteio*raio*cos(angulo);
	    // p[i].y=L/2.+raiosorteio*raio*sin(angulo);
	    p[i].x=L*ran2();
	    p[i].y=L*ran2();
	     // sorteio inicial das velocidades	
	     angle=TWOPI*ran2();
	     p[i].theta=angle;
	     p[i].vx=p[i].v0*cos(angle);
	     p[i].vy=p[i].v0*sin(angle);
	  }
	cudaMemcpy(Gp,p,memP,cudaMemcpyHostToDevice);  
     }   
    //since we need the particle size and force radius only in squared form we do:
    a*=a;
    ra*=ra;
    return 0;
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	particles *pc=new particles();
	pc->initparams(argc,argv);
	printf("params set, initializing ... ");
	if(pc->initialize()!=0)
	{
	   printf("program stopped\n");
	   return 0;
	}
	printf("done\nstarting the simulation\n");
	if(pc->simtype==0) pc->simulate();
	printf("finished - cleaning up\n");
	delete pc;
	return 0;
}
//---------------------------------------------------------------------------



