#include <cuda.h>
#include <stdlib.h>
using namespace std;

extern "C" {
  void pninpoly(int nvert, float *vertx, float *verty, int ntest, \
		float *testx, float *testy, int *mask);
} 

//polong inside of a polygon
__global__ void pnpoly(const int nvert, const float * __restrict__ vertx, const float * __restrict__ verty, const int ntest, \
		       const float * __restrict__ testx, const float * __restrict__ testy, int *mask){
  int i,j,c;
  int tid;
  tid=threadIdx.x+blockIdx.x*blockDim.x;
  while(tid<ntest){
    c=0;
    for(i=0,j=nvert-1;i<nvert;j=i++){
      if((((verty[i]<=testy[tid]) &&(testy[tid]<verty[j]))||
	  ((verty[j]<=testy[tid]) &&(testy[tid]<verty[i])))&&
	 (testx[tid]<(vertx[j]-vertx[i])*(testy[tid]-verty[i])/(verty[j]-verty[i])+vertx[i])) c=!c;
    }
    mask[tid]=c;
    tid+=blockDim.x*gridDim.x;
  }
}
//
void pninpoly(int nvert, float *vertx, float *verty, int ntest, \
	      float *testx, float *testy, int *mask){
  float *d_vx, *d_vy, *d_tx, *d_ty;
  int *d_mk;
  
  cudaMalloc( (void**)&d_vx, nvert*sizeof(float));
  cudaMalloc( (void**)&d_vy, nvert*sizeof(float));
  cudaMalloc( (void**)&d_tx, ntest*sizeof(float));
  cudaMalloc( (void**)&d_ty, ntest*sizeof(float));
  cudaMalloc( (void**)&d_mk, ntest*sizeof(int));
  //cudaMalloc( (void**)&d_dy,sizeof(float));

  //copy to gpu
  cudaMemcpy(d_vx,vertx,nvert*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_vy,verty,nvert*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_tx,testx,ntest*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_ty,testy,ntest*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_mk,mask,ntest*sizeof(int),cudaMemcpyHostToDevice);
  //cudaMemcpy(d_dy,dy,sizeof(float),cudaMemcpyHostToDevice);

  pnpoly<<<1024,1024>>>(nvert,d_vx,d_vy,ntest,d_tx,d_ty,d_mk);

  cudaMemcpy(mask,d_mk,ntest*sizeof(int),cudaMemcpyDeviceToHost);
  cudaFree(d_vx);cudaFree(d_vy);cudaFree(d_tx);cudaFree(d_ty);cudaFree(d_mk);
}

  
