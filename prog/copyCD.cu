#include <cmath>
#include <cuda.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

extern "C" {
  void copyCDs(int nx, float *cdx, int ny, float *cdy, int nz, float *cdz, float * cdose2d, \
		       int npcd3, float *rax, float *ray, float *raz, float * cDose3d,
		       float dx, float dy, float dz);
} 


__device__ int wherb(const float * __restrict__ x,const float val,const int n){
  int i,indx=0;
  float minn;
  minn=(abs(x[0]-val));

  for(i=1;i<n;i++){
    if (abs(x[i]-val)< minn){
      minn=abs(x[i]-val);
      indx=i;
    }
  }
  return indx;
}

__device__ int wher(const float *__restrict__ x,const float val, const int n){
  int nabove=n+1,nbelow=0,mid;
  //pass in dlo/2
    while(nabove -nbelow >1){
      mid=(nabove+nbelow)/2;
      if(val==x[mid-1])return mid-1;
      if(val< x[mid-1]){ nabove=mid;}
      else{ nbelow=mid;}
    }
    return nbelow-1;
  }

//trilinear interpolation
//
__global__ void copyCD(const int nx, const float *__restrict__ cdx, const int ny, const float * __restrict__ cdy, \
                const int nz, const float * __restrict__ cdz, float * cdose2d, \
		       const int npcd3, const float * __restrict__ rax, const float * __restrict__ ray, const float * __restrict__ raz,\
                const float * __restrict__ cDose3d,
		       const float dx, const float dy, const float dz){
  int i,j,k;
  int tid;
  tid=threadIdx.x+blockIdx.x*blockDim.x;
  while(tid<npcd3){
    i=wherb(cdy,ray[tid],ny);
    j=wherb(cdx,rax[tid],nx);
    k=wherb(cdz,raz[tid],nz);
    if(abs(cdy[i]-ray[tid])<=dy/2.0 && abs(cdx[j]-rax[tid])<=dx/2.0 || abs(cdz[k]-raz[tid])<=dz/2.0){
      cdose2d[k+j*nz+i*nz*nx]=cDose3d[tid];
    }else{
      cdose2d[k+j*nz+i*nz*nx]=0.0f;
    }
    tid+=blockDim.x*gridDim.x;
  }
}
//
void copyCDs(int nx, float *cdx, int ny, float *cdy, int nz, float *cdz, float * cdose2d, \
		       int npcd3, float *rax, float *ray, float *raz, float * cDose3d,
		       float dx, float dy, float dz){
  float *d_cdx, *d_cdy, *d_cdz, *d_cdos2;
  float *d_rax, *d_ray, *d_raz, *d_cdos3;
  //float *h_cdos2,*h_cdos3;
  //float *h_rax, *h_ray, *h_raz, *h_cdos3;
  //int i,curi;
  cudaMalloc( (void**)&d_cdx, nx*sizeof(float));
  cudaMalloc( (void**)&d_cdy, ny*sizeof(float));
  cudaMalloc( (void**)&d_cdz, nz*sizeof(float));
  cudaMalloc( (void**)&d_cdos2,nz*nx*ny*sizeof(float));
  cudaMalloc( (void**)&d_rax, npcd3*sizeof(float));
  cudaMalloc( (void**)&d_ray, npcd3*sizeof(float));
  cudaMalloc( (void**)&d_raz, npcd3*sizeof(float));
  cudaMalloc( (void**)&d_cdos3,npcd3*sizeof(float));
  //cudaHostAlloc( (void**)&h_cdx, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_cdy, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_cdz, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_cdos2,npcd2*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_cdos3,npcd3*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_rax, rx*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_ray, ry*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_raz, rz*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //cudaHostAlloc( (void**)&h_cdos3,rz*rx*ry*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  //for(i=0;i<npcd2;i++){
  //  h_cdos2[i]=cdose2d[i];
  // }
  //for(i=0;i<npcd3;i++){
  //  h_cdos3[i]=cDose3d[i];
  //}
  
  //copy to gpu
  cout<<cudaMemcpy(d_cdx,cdx,nx*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_cdy,cdy,ny*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_cdz,cdz,nz*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_cdos2,cdose2d,nz*nx*ny*sizeof(float),cudaMemcpyHostToDevice)<<endl;;
  cout<<cudaMemcpy(d_rax,rax,npcd3*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_ray,ray,npcd3*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_raz,raz,npcd3*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  cout<<cudaMemcpy(d_cdos3,cDose3d,npcd3*sizeof(float),cudaMemcpyHostToDevice)<<endl;
  //cout<<cudaHostGetDevicePointer(&d_cdos2,h_cdos2,0)<<endl;
  //cout<<cudaHostGetDevicePointer(&d_cdos3,h_cdos3,0)<<endl;
  //curi=48*128;
  //for(i=0;i<npcd2;i+=curi){
  //  cout<<"percent done "<<float(i)/float(npcd2)*100<<endl;
    
  copyCD<<<512,512>>>(nx,d_cdx,ny,d_cdy,nz,d_cdz,d_cdos2,	\
		       npcd3,d_rax,d_ray,d_raz,d_cdos3,	\
		       dx,dy,dz);
  //}

  cout<<cudaMemcpy(cdose2d,d_cdos2,nz*nx*ny*sizeof(float),cudaMemcpyDeviceToHost)<<endl;
  //for(i=0;i<npcd2;i++){
  //  cdose2d[i]=h_cdos2[i];
  //}
  cout<<"finished"<<endl;
  //cudaDeviceReset();
  cudaFree(d_cdx);cudaFree(d_cdy);cudaFree(d_cdz);cudaFree(d_cdos2);
  cudaFree(d_rax);cudaFree(d_ray);cudaFree(d_raz);cudaFree(d_cdos3);
  //cudaFreeHost(h_cdos2);cudaFreeHost(h_cdos3);
}

