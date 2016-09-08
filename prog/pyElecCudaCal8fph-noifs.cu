#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cuda.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
//#include <thrust/device_ptr.h>
//#include <thrust/device_delete.h>
//#include <thrust/sort.h>
//#include <thrust/extrema.h>
using namespace std;
/*Used to interact with python*/
extern "C" {
  void pymain (float E0,int NROW,int NCOL,int NZ, int npcell, int nMP, int nthreads, float alp, \
	       float A, float bnc, int dens[],	\
	       float x[], float y[], float z[], float dose[], float strt, \
	       float sstop, float xmin, float xmax,float blkx[], float blky[], \
	       int nblk, float ymin, float ymax,\
	       float zmin, float zmax,float phper, float felec[],\
	       float sx, float ix, float sy, float iy, float sz, float iz);
}
//global variables
//const int NMAX=75;
const float PI=3.1415926535897932385f;

#include "incs2f-noifs.h"
//sorting functions
__global__ void bitonic_sort_step(float *dev_values, int j, int k)
{
  unsigned int i, ixj; /* Sorting partners: i and ixj */
  i = threadIdx.x + blockDim.x * blockIdx.x;
  ixj = i^j;

  /* The threads with the lowest ids sort the array. */
  if ((ixj)>i) {
    if ((i&k)==0) {
      /* Sort ascending */
      if (dev_values[i]>dev_values[ixj]) {
        /* exchange(i,ixj); */
        float temp = dev_values[i];
        dev_values[i] = dev_values[ixj];
        dev_values[ixj] = temp;
      }
    }
    if ((i&k)!=0) {
      /* Sort descending */
      if (dev_values[i]<dev_values[ixj]) {
        /* exchange(i,ixj); */
        float temp = dev_values[i];
        dev_values[i] = dev_values[ixj];
        dev_values[ixj] = temp;
      }
    }
  }
}

/**
 * Inplace bitonic sort using CUDA.
 */
void bitonic_sort(float *values,const int nMP, const int nthreads)
{
  //float *dev_values;
  size_t size = nMP*nthreads * sizeof(float);

  //cudaMalloc((void**) &dev_values, size);
  //cudaMemcpy(dev_values, values, size, cudaMemcpyHostToDevice);

  dim3 blocks(nMP,1);    /* Number of blocks   */
  dim3 threads(nthreads,1);  /* Number of threads  */

  int j, k;
  /* Major step */
  for (k = 2; k <= nMP*nthreads; k <<= 1) {
    /* Minor step */
    for (j=k>>1; j>0; j=j>>1) {
      bitonic_sort_step<<<blocks, threads>>>(values, j, k);
    }
  }
  //cudaMemcpy(values, dev_values, size, cudaMemcpyDeviceToHost);
  //cudaFree(dev_values);
}

__device__ float delSr(const float cfac,const float sfr,const float sfc, const float e0, const float tdl){
	/*change in energy by moving through with radiative interaction SeSr*/
	/*CUDA for delE=fac*sfr*SeSr(e0)*tdl;*/
	return __fmul_rn(__fmul_rn(__fmul_rn(cfac,sfr),SeSr(e0)),tdl);
}
__device__ float delSc(const float cfac,const float sfr,const float sfc, const float e0, const float tdl){
	/*change in energy by moving through with collisional interaction SeSc*/
	/*CUDA for delE=sfc*SeSc(e0)*tdl*fac;*/
	return __fmul_rn(__fmul_rn(__fmul_rn(sfc,SeSc(e0)),tdl),cfac);
}

//my rng
__device__ inline float xor128(uint *state)
{
state[0]=36969*(state[0]&65535)+(state[0]>>16)+69069*state[0]+1234567;
state[1]=18000*(state[1]&65535)+(state[1]>>16);
return ((state[0]<<16)+state[1])*2.328306e-10f;
}

__device__ inline float normR(const float u1, const float u2){
        return sqrt(-2.0f*__logf(u1))*__cosf(2.0f*PI*u2);
}

__global__ void normE0s( float *e0s, uint *state0, uint *state1, const float alp, const float E0){
  int init=threadIdx.x+blockIdx.x*blockDim.x;
  uint localstate[2];
  if(init<gridDim.x*blockDim.x){
    localstate[0]=state0[init];localstate[1]=state1[init];
    /*get 2 rngs*/
    //float r1=xor128(localstate);
    //float r2=xor128(localstate);
    /*get random normal*/
    e0s[init]=normR(xor128(localstate),xor128(localstate))*alp+E0;
    /*update global rng state*/
    state0[init]=localstate[0];state1[init]=localstate[1];
  }
}
__global__ void pos( float *exs, float *ezs, uint *state0, uint *state1, const float sstp, const float strt){
  int init=threadIdx.x+blockIdx.x*blockDim.x;
  uint localstate[2];
  if(init<gridDim.x*blockDim.x){
    localstate[0]=state0[init];localstate[1]=state1[init];
    /*generate uniform numbers from strt to stop*/
    exs[init]=xor128(localstate)*(sstp-strt)-(sstp-strt)/2.0f;
    //exs[init]=(sstp-strt)*threadIdx.x/blockDim.x+strt;
    ezs[init]=xor128(localstate)*(sstp-strt)-(sstp-strt)/2.0f;
    /*update global rng state*/
    state0[init]=localstate[0];state1[init]=localstate[1];
  }
}
__global__ void npnpoly( int * ntins, const int nvert,const float * __restrict__ vertx, \
			  const float * __restrict__ verty, const float * __restrict exs, \
			  const float * __restrict__ ezs){
  int init=threadIdx.x+blockIdx.x*blockDim.x;
  if(init <gridDim.x*blockDim.x){
    ntins[init]=pnpoly(nvert,vertx,verty,exs[init],ezs[init]);
  }
}
/*Main calculation engine*/
__global__ void kernel(  const float dxyz,float * dose,			\
			 const int * __restrict__ dens,			\
			 uint *state0, uint *state1,const float * __restrict__ e0s, const float * __restrict__ exs, const float * __restrict__ ezs, \
			 const int * __restrict__ ntins, const float A,	\
			 const float bnc,				\
			 const int NROW, const int NCOL, const int NZ,	\
			 const float xmin, const float xmax, const float ymin, \
			 const float ymax,const float zmin, const float zmax,const float phper, \
			 const float sx, const float ix, const float sy, const float iy, const float sz, const float iz){
  //const float * __restrict__ x, const float * __restrict__ y, const float * __restrict__ z,
  float minSe=0.01f,e0,delE;
  float ex, ey, ez, edx, edy, edz,vssd=90.0f;
  float theta,the1,the2,urn,dlo,r1,fac,afac=1.2e-3f,cfac,sfc,sfr,sfy,fx,fy,fz;
  float ddens[9];
  int idens;
  ddens[0]=1.20e-3f;//air
  ddens[1]=0.2f;ddens[2]=0.5f;//ilung, elung
  ddens[3]=0.96f;ddens[4]=1.0f;//fat, water
  ddens[5]=1.04f;ddens[6]=1.15f;//muscle, cbone
  ddens[7]=1.85f;ddens[8]=1.85f;//sbone, lead (but bone)
  float vec0,vec1,vec2,tdl=0.0f;//,evecx[2],evecy[2],evecz[2];
  uint localstate[2];
  int indx=0,indy=0,indz=0,notin;
  int init,ph=0;
  /*set electron translation delta, bnc forces multiple scattering in the same voxel*/
  dlo=dxyz*bnc;
  /*defines the unique thread number using CUDA variables*/
  init=threadIdx.x+blockIdx.x*blockDim.x;
  /*choose the local random number generator state*/
  localstate[0]=state0[init];localstate[1]=state1[init];
  /*get electron initial values*/
  /*initial energy*/
  e0=e0s[init];
  
  /*initial position in x, presorted to help with divergent executions*/
  ex=exs[init];
  /*initial position in z*/
  ez=ezs[init];
  /*initial position in y*/
  ey=ymin+dlo;
  /*is the elctron in the grid and in the apeture of the block, precomputed before kernel call*/
  notin=ntins[init];
  //if(init==0){
  //  e0=12.0f;ex=0.0f;ez=0.0f;notin=1.0f;
  //}
  //if (1-notin) e0=0.0f;
  e0*=(float)notin;
  /*magnitude of {ex,ez,vssd} used to scale the step sizes {edx,edy,edz}*/
  urn=sqrt(ex*ex+ez*ez+vssd*vssd);
  /*stepsize pseudo velocity vector, used to step the electrons through the grid*/
  edy=dlo*vssd/urn;edx=dlo*ex/urn;edz=dlo*ez/urn;
  //is this a photon
  /*phper from python is the photon percentage needed to match the machine measurements is random number is less than the percentage then make this a photon and not an electron*/
  //if(xor128(localstate)<phper)ph=1;
  ph=signbit(xor128(localstate)-phper);
  /*begin energy deposition*/
  //if(ez>zmax || ez<zmin||ex>xmax||ex<xmin||ey>ymax||ey<ymin)e0=0.0f;
  while(e0>minSe){
    /*if the electron has wandered off the grid exit the loop*/
    //checked below indirectly
    //if(ez>zmax || ez<zmin||ex>xmax||ex<xmin||ey>ymax||ey<ymin)break;
    /*find which array index the electron is in using binary search, only search if it has moved out of the current voxel, see include file*/
    //indx=d_wherb(x,ex,NCOL,indx);indy=d_wherb(y,ey,NROW,indy);
    //indz=d_wherb(z,ez,NZ,indz);      
    indx=(int)floor((ex-ix)/sx);indy=(int)floor((ey-iy)/sy);indz=(int)floor((ez-iz)/sz);
    /*fx=(float)signbit(xmax-ex)+(float)signbit(ex-xmin);//+signbit((float)(NCOL-1)-indx));
    fy=(float)signbit(ymax-ey)+(float)signbit(ey-ymin);//+signbit((float)(NROW-1)-indy));
    fz=(float)signbit(zmax-ez)+(float)signbit(ez-zmin);//+signbit((float)(NZ-1)-indz));
    e0=(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*e0;*/
    /*fx=(float)signbit(ex-xmin);
    fy=(float)signbit(ey-ymin);
    fz=(float)signbit(ez-zmin);
    e0=(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*e0;*/
    //if(init==0)printf("indx %d, indy %d, indz %d\n",indx,indy,indz);
    //if(indx < 0 || indx > NCOL-1 || indy < 0 || indy > NROW-1 || indz < 0 || indz > NZ-1)break;
    //if(init==0)printf("indx %d, indy %d, indz %d\n",indx,indy,indz);
    //check the density we are currently in an integer array sent by python 
    idens=dens[indz+indx*NZ+indy*NZ*NCOL];
    fac=ddens[idens];
    /*scale the water stopping powers etc to the current material*/
    sfc=scalSeSc(e0,idens);
    sfr=scalSeSr(e0,idens);
    sfy=scalSeSy(e0,idens);
    /*scatter the electron for next iteration*/
    vec0=edx;vec1=edy;vec2=edz;
    /*adjust for photon vs electron*/
    //if(ph){cfac=afac;}else{cfac=fac;}
    cfac=__fmul_rn(afac,(float)ph)+__fmul_rn(fac,(1.0f-(float)ph));
    /*CUDA notation for theta=sqrt(fac*dl)*(A*PI/(180))*/
    theta=__fmul_rn(__fsqrt_rn(__fmul_rn(cfac,dlo)),__fdividef(__fmul_rn(A,PI),180.0f));
    /*for photon or air density scattring angle is larger by 1.5 times, experimentally determined*/
    //if(fac==afac){theta*=0.5f;}
    /* random number to scale the scattering angle from (-theta,theta)*/
    urn=__fadd_rn(__fmul_rn(xor128(localstate),2.0f),-1.0f);
    /*find the direction with the largest magnitude, scater around it randomly choosing one of the other two directions*/
    theta=__fmul_rn(theta,urn);
    the1=__cosf(theta);the2=__sinf(theta);
    /*random number for decision tree below*/
    float fdesc=round(xor128(localstate));
    float omfdesc=1.0f-fdesc;
    //sedx=fabs(edx);sedy=fabs(edy);sedz=fabs(edz);
    fx=(float)signbit(edy-edx)*signbit(edz-edx);
    fy=(float)signbit(edx-edy)*signbit(edz-edy);
    fz=(float)signbit(edx-edz)*signbit(edy-edz);
    //if(edy >= edx && edy>=edz){
      //if(idesc){
    edx=(vec0*fdesc+(__fadd_rn(__fmul_rn(vec0,the1),-__fmul_rn(vec1,the2)))*omfdesc)*fy+\
      ((__fadd_rn(__fmul_rn(vec0,the1),__fmul_rn(vec2,the2)))*fdesc+(__fadd_rn(__fmul_rn(vec0,the1),-__fmul_rn(vec1,the2)))*omfdesc)*fx+\
      (vec0*fdesc+(__fadd_rn(__fmul_rn(vec0,the1),__fmul_rn(vec2,the2)))*omfdesc)*fz;
    edy=((__fadd_rn(__fmul_rn(vec1,the1),-__fmul_rn(vec2,the2)))*fdesc+(__fadd_rn(__fmul_rn(vec1,the1),__fmul_rn(vec0,the2)))*omfdesc)*fy+\
      (vec1*fdesc+(__fadd_rn(__fmul_rn(vec1,the1),__fmul_rn(vec0,the2)))*omfdesc)*fx+\
      ((__fadd_rn(__fmul_rn(vec1,the1),-__fmul_rn(vec2,the2)))*fdesc+vec1*omfdesc)*fz;
    edz=((__fadd_rn(__fmul_rn(vec2,the1),__fmul_rn(vec1,the2)))*fdesc+vec2*omfdesc)*fy+\
      ((__fadd_rn(__fmul_rn(vec2,the1),-__fmul_rn(vec0,the2)))*fdesc+vec2*omfdesc)*fx+\
      ((__fadd_rn(__fmul_rn(vec2,the1),__fmul_rn(vec1,the2)))*fdesc+(__fadd_rn(__fmul_rn(vec2,the1),-__fmul_rn(vec0,the2)))*omfdesc)*fz;
	//}else{
	//edx=__fadd_rn(__fmul_rn(vec0,the1),-__fmul_rn(vec1,the2));
	//edy=__fadd_rn(__fmul_rn(vec1,the1),__fmul_rn(vec0,the2));
      //}
      //}else if(edx >= edy && edx >= edz){
    //if(idesc){
	/*CUDA for 
	  edx=(vec[0]*cos(the2)+vec[2]*sin(the2));
	  edz=(vec[2]*cos(the2)-vec[0]*sin(the2));*/
    //edx+=((__fadd_rn(__fmul_rn(vec0,the1),__fmul_rn(vec2,the2)))*fdesc+(__fadd_rn(__fmul_rn(vec0,the1),-__fmul_rn(vec1,the2)))*omfdesc)*fx;
    //edy+=(vec1*fdesc+(__fadd_rn(__fmul_rn(vec1,the1),__fmul_rn(vec0,the2)))*omfdesc)*fx;
    //edz+=((__fadd_rn(__fmul_rn(vec2,the1),-__fmul_rn(vec0,the2)))*fdesc+vec2*omfdesc)*fx;
      //}else{
	/*CUDA for 
	  edx=(vec[0]*cos(the3)-vec[1]*sin(the3));
	  edy=(vec[1]*cos(the3)+vec[0]*sin(the3));*/
	//edx=__fadd_rn(__fmul_rn(vec0,the1),-__fmul_rn(vec1,the2));
	//edy=__fadd_rn(__fmul_rn(vec1,the1),__fmul_rn(vec0,the2));
      //}
      /*the above pattern continues below for edy and edz*/
      //}else if(edz>=edx && edz>=edy){
      //if(idesc){
    //edx+=(vec0*fdesc+(__fadd_rn(__fmul_rn(vec0,the1),__fmul_rn(vec2,the2)))*omfdesc)*fz;
    //edy+=((__fadd_rn(__fmul_rn(vec1,the1),-__fmul_rn(vec2,the2)))*fdesc+vec1*omfdesc)*fz;
    //edz+=((__fadd_rn(__fmul_rn(vec2,the1),__fmul_rn(vec1,the2)))*fdesc+(__fadd_rn(__fmul_rn(vec2,the1),-__fmul_rn(vec0,the2)))*omfdesc)*fz;
	//}else{
	//edx=__fadd_rn(__fmul_rn(vec0,the1),__fmul_rn(vec2,the2));
	//edz=__fadd_rn(__fmul_rn(vec2,the1),-__fmul_rn(vec0,the2));
	//}
    //}
    /*move the elctron to the next position*/
    ex+=(edx+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec0);    
    ey+=(edy+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec1);
    ez+=(edz+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec2);

    //prevent 0 being propagated
    edx=(edx+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec0);
    edy=(edy+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec1);
    edz=(edz+(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*vec2);
    /*if(init==0)printf("edx= %g, edy= %g, edz= %g\n",edx,edy,edz);
      if(init==0)printf("vec0= %g, vec1= %g, vec2= %g\n",vec0,vec1,vec2);*/
    //if(init==0)printf("fx= %g, fy= %g, fz= %g\n",ex,ey,ez);
    /*total distance moved, CUDA for sqrt(edx*edx+edy*edy+edz*edz);*/
    tdl=__fsqrt_rn(__fmul_rn(edx,edx)+__fmul_rn(edy,edy)+__fmul_rn(edz,edz));
    /*CUDA forr1=sfy*SeSy(e0); choose wheter to have a radiative interaction of collision interaction*/ 
    r1=__fmul_rn(sfy,SeSy(e0));
    /*fx is used to get rid of the if statement*/
    fx=(float)signbit(xor128(localstate)-r1);
    /*how much energy is lost*/
    delE=delSc(cfac,sfr,sfc,e0,tdl)*(1.0f-fx)+delSr(cfac,sfr,sfc,e0,tdl)*fx;
    //if(init==0)printf("fx= %f dele= %f\n",fx,delE);
    fx=(float)signbit(e0-delE);
    e0-=(1.0f-fx)*delE+fx*e0;
    //if(delE>e0)delE=e0;
    //e0-=delE;
    /*add the lost energy to the photon dose array*/
    //if(init==0)printf("e0 = %f urn= %f tdl= %f sfc= %f SeSc= %f sfr= %f SeSr=%f\n",e0,delE/cfac,tdl,sfc,SeSc(e0),sfr,SeSr(e0));
    atomicAdd(&dose[indz+indx*NZ+indy*NZ*NCOL],__fdividef(delE,cfac));
    fx=(float)signbit(xmax-ex)+(float)signbit(ex-xmin);//+signbit((float)(NCOL-1)-indx));
    fy=(float)signbit(ymax-ey)+(float)signbit(ey-ymin);//+signbit((float)(NROW-1)-indy));
    fz=(float)signbit(zmax-ez)+(float)signbit(ez-zmin);//+signbit((float)(NZ-1)-indz));
    e0=(1.0f-fx)*(1.0f-fy)*(1.0f-fz)*e0;
    //indx=(int)floor((0.0f-ix)/sx);indy=(int)floor((2.85f-iy)/sy);indz=(int)floor((0.0f-iz)/sz);
    //if(init==0)printf("indx= %d indy= %d indz=%d dose= %f\n",indx,indy,indz,dose[indz+indx*NZ+indy*NZ*NCOL]);
  }/*ends while loop*/
  /*save the local random number state for the next call to the kernel*/
  state0[init]=localstate[0];state1[init]=localstate[1];
}
  
void pymain (float E0,int NROW, int NCOL, int NZ, int npcell, int nMP, int nthreads, float alp, \
	     float A, float bnc, int dens[],				\
	     float x[], float y[], float z[], float dose[], float strt, \
	     float sstop, float xmin, float xmax,float blkx[], float blky[], \
	     int nblk, float ymin, float ymax,				\
	     float zmin, float zmax,float phper,float felec[],\
	     float sx, float ix, float sy, float iy, float sz, float iz){
  //variable for random numbers
  const gsl_rng_type * T;
  gsl_rng * r;
  /*let compiler know none of these will change*/
  const float AA=A,bbnc=bnc,EE0=E0,aalp=alp,sstrt=strt,sstp=sstop,xxmin=xmin,xxmax=xmax;
  const float yymin=ymin,yymax=ymax,zzmin=zmin,zzmax=zmax,pphper=phper,ssx=sx,iix=ix,ssy=sy,iiy=iy,ssz=sz,iiz=iz;
  const int NNROW=(int)NROW,NNCOL=(int)NCOL,NNZ=(int)NZ,nnblk=(int)nblk,nnMP=(int)nMP,nnthreads=(int)nthreads,nnpcell=(int)npcell;
  const float dxyz=fabs(x[1]-x[0]);
  /*pointers to the CUDA device d_ and host h_ arrays*/
  //float *d_x,*d_y,*d_z;
  float *d_blkx, *d_blky;
  float *d_dose,*d_e0s, *d_exs, *d_ezs;
  int *d_dens, *d_ntins;
  uint *d_seeds0,*d_seeds1;
  /*pointers for host arrays*/
  uint *seeds0,*seeds1;
  float *tdose;
  int i,k;
  int indx1,indx2,indy1,indy2;
  int totelec=0;
  int nelec,nloop;
  cout<<"E0= "<<E0<<" A= "<<A<<" alp= "<<alp<<" bnc= "<<bnc<<endl;
  cout<<"NROW= "<<NROW<<" NCOL= "<<NCOL<<" NZ= "<<NZ<<endl;
  cout<<"Begin PyMain"<<endl;
  //init random number generator
  gsl_rng_env_setup();
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
  gsl_rng_set(r,random_seed());
  //setup domain
  tdose= new float[NNROW*NNCOL*NNZ];  
  /*initialize arrays*/
  for(i=0;i<NNROW*NNCOL*NNZ;i++){
    tdose[i]=0.0f;
  }
  //initalize random seeds
  seeds0=new uint[nnMP*nnthreads];
  seeds1=new uint[nnMP*nnthreads];
  /*count the number of pixels in a beams eye view and make the total electron count
    number of pixels * electrons per pixel from python*/
  indx1=wherb(x,sstrt,NNCOL);indx2=wherb(x,sstp,NNCOL);
  indy1=wherb(z,sstrt,NNZ);indy2=wherb(z,sstp,NNZ);
  /*get the number of electrons*/
  nelec=int((indx2-indx1+1))*int((indy2-indy1+1))*int(nnpcell);
  //redefine NPCELL to loop over kernel parameter
  nloop=(nelec/(nnMP*nnthreads))+1;
  cout<<"nelec = "<<nelec<<" width= "<<(indx2-indx1+1)<<" height= "<<(indy2-indy1+1)<<endl;
  cout<<"nloops= "<<nloop<<endl;
  /*allocate cuda arrays*/
  //cout<<"alloc x"<<endl;
  /*allocate the device memory for x,y,z,dose,dens,blocks*/
  //cudaMalloc( (void**)&d_x,NCOL*sizeof(float));
  /*copy from host to device*/
  //cudaMemcpy(d_x,x,NCOL*sizeof(float),cudaMemcpyHostToDevice);
  //cout<<"alloc y"<<endl;
  //cudaMalloc( (void**)&d_y,NROW*size(float));
  //cudaMemcpy(d_y,y,NROW*sizeof(float),cudaMemcpyHostToDevice);
  //cout<<"alloc z"<<endl;
  //cudaMalloc( (void**)&d_z,NZ*sizeof(float));
  //cudaMemcpy(d_z,z,NZ*sizeof(float),cudaMemcpyHostToDevice);
  cout<<"alloc dens"<<endl;
  cudaMalloc( (void**)&d_dens,NNROW*NNCOL*NNZ*sizeof(int));
  cudaMemcpy(d_dens,dens,NNROW*NNCOL*NNZ*sizeof(int),cudaMemcpyHostToDevice);
  cout<<"alloc dose"<<endl;
  cudaMalloc( (void**)&d_dose,NNROW*NNCOL*NNZ*sizeof(float));
  cudaMemcpy(d_dose,dose,NNROW*NNCOL*NNZ*sizeof(float),cudaMemcpyHostToDevice);
  cout<<"alloc e0s"<<endl;
  cudaMalloc( (void**)&d_e0s,nnMP*nnthreads*sizeof(float));
  cudaMalloc( (void**)&d_exs,nnMP*nnthreads*sizeof(float));
  cudaMalloc( (void**)&d_ezs,nnMP*nnthreads*sizeof(float));
  cudaMalloc( (void**)&d_ntins,nnMP*nnthreads*sizeof(int));
  /*thrust pointer to exs*/
  //thrust::device_ptr<float> dev_ptr(d_exs);
  //thrust::device_ptr<float> dev_ptr0(d_e0s);
  
  cout<<"alloc blkx"<<endl;
  cudaMalloc((void**)&d_blkx,nnblk*sizeof(float));
  cudaMemcpy(d_blkx,blkx,nnblk*sizeof(float),cudaMemcpyHostToDevice);
  cout<<"alloc blky"<<endl;
  cudaMalloc((void**)&d_blky,nnblk*sizeof(float));
  cudaMemcpy(d_blky,blky,nnblk*sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_seeds0,nnMP*nnthreads*sizeof(uint));
  cudaMalloc((void**)&d_seeds1,nnMP*nnthreads*sizeof(uint));
  for(i=0;i<nnMP*nnthreads;i++){
    seeds0[i]=gsl_rng_get(r);
    seeds1[i]=gsl_rng_get(r);
  }
  cudaMemcpy(d_seeds0,seeds0,nnMP*nnthreads*sizeof(uint),cudaMemcpyHostToDevice);
  cudaMemcpy(d_seeds1,seeds1,nnMP*nnthreads*sizeof(uint),cudaMemcpyHostToDevice);
  //redefine nelec here WARNING
  nelec=nnMP*nnthreads;
  cout<<"nelec redefine= "<<nelec<<endl;
  //thrust::device_ptr<float> max_ptr=thrust::max_element(dev_ptr0,dev_ptr+nelec);
  for (k=0;k<nloop;k++){
    cout<<"Kernel Loop "<<float(k)*100.0/float(nloop)<<" % done"<<endl;
    /*call the kernel until all electrons have been calculated*/
    /*generate random normal distribution for electron energies*/
    normE0s<<<nnMP,nnthreads>>>(d_e0s,d_seeds0,d_seeds1,aalp,EE0);
    /*initialize electron positions*/
    pos<<<nnMP,nnthreads>>>(d_exs,d_ezs,d_seeds0,d_seeds1,sstop,sstrt);
    /*sort exs with thrust*/
    //thrust::sort(dev_ptr,dev_ptr+nelec);
    bitonic_sort(d_exs,nnMP,nnthreads);
    /*are they in the block*/
    npnpoly<<<nnMP,nnthreads>>>(d_ntins,nnblk,d_blkx,d_blky,d_exs,d_ezs);
    /*calculate dose using precomputed and sorted electrons*/
    //while(max_ptr[0]>0.01f){
    kernel<<<nnMP,nnthreads>>>(dxyz,d_dose,d_dens,			\
			     d_seeds0,d_seeds1,d_e0s,d_exs,d_ezs,d_ntins,AA,bbnc,NNROW,NNCOL,NNZ, \
			     xxmin,xxmax,yymin,yymax,zzmin,zzmax,pphper,ssx,iix,ssy,iiy,ssz,iiz);
    //  max_ptr=thrust::max_element(dev_ptr0,dev_ptr+nelec);
    //}
    cudaThreadSynchronize();
    /*what is happening with dose*/
    //cudaMemcpy(tdose,d_dose,NNROW*NNCOL*NNZ*sizeof(float),cudaMemcpyDeviceToHost);
    //cout<<"MAX dose= "<<maxval(tdose,NNROW*NNCOL*NNZ)<<endl;
    /*count total actual electrons*/
    totelec+=nelec;
  }
  /*copy doses back from device*/
  cudaMemcpy(tdose,d_dose,NNROW*NNCOL*NNZ*sizeof(float),cudaMemcpyDeviceToHost);
  cout<<"nelec tot= "<<totelec<<endl;
  /*update python dose array*/
  for(i=0;i<NNROW*NNCOL*NNZ;i++){
    dose[i]=tdose[i];//*1.602e-10;
  }
  felec[0]=static_cast<float>(totelec);
  cout<<"MAX dose= "<<maxval(dose,NNROW*NNCOL*NNZ)<<endl;
  gsl_rng_free(r);
  delete seeds0; delete seeds1;
  delete tdose;
  //thrust::device_delete(dev_ptr);
  //thrust::device_delete(dev_ptr0);
  cudaFree(d_dens);cudaFree(d_dose);cudaFree(d_e0s);cudaFree(d_exs);cudaFree(d_ezs);
  cudaFree(d_ntins);cudaFree(d_blkx);cudaFree(d_blky);cudaFree(d_seeds0);cudaFree(d_seeds1);
  //cudaFree(dev_ptr);
  cudaDeviceReset();
}
