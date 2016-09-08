#include <cmath>
#include <cuda.h>
#include <stdlib.h>

using namespace std;

extern "C" {
  void interps(int npnts, float *caax, float *caay, float *caaz,float *cDose3dr, \
	       int NX, float *cax, int NY, float *cay,			\
	       int NZ, float *caz, float * cDose3d,float * mnmxs);
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
__global__ void interp(const int npnts, const float *__restrict__ caax, const float * __restrict__ caay, const float * __restrict__ caaz,float *cDose3dr, \
		       const int NX, const float * __restrict__ cax, const int NY,const float * __restrict__ cay,\
		       const int NZ, const float * __restrict__ caz, const float * __restrict__ cDose3d,const float * __restrict__ mnmxs){
  //mnmxs passes in minimum and maximum of each array y->0,1 x->2,3 z->4,5
  int i0,i1,j0,j1,k0,k1;
  int tid=threadIdx.x+blockIdx.x*blockDim.x;;
  float xc,yc,zc;
  float c00,c10,c01,c11,c0,c1;
  float cx1,cx0,cy0,cy1,cz0,cz1;
  while(tid<npnts){
    i0=wherb(cay,caay[tid],NY);
    if(i0==NY-1)i0=i0-1;
    i1=i0+1;
    j0=wherb(cax,caax[tid],NX);
    if(j0==NX-1)j0=j0-1;
    j1=j0+1;
    k0=wherb(caz,caaz[tid],NZ);
    if(k0==NZ-1)k0=k0-1;
    k1=k0+1;
    if(caay[tid]<mnmxs[0] || caay[tid] >mnmxs[1] || caax[tid]<mnmxs[2] ||	\
       caax[tid]> mnmxs[3] || caaz[tid]<mnmxs[4] || caaz[tid]>mnmxs[5] || \
       cay[i1]==cay[i0] || cax[j1]==cax[j0] || caz[k1]==caz[k0] || \
       i0<0 || j0<0 || k0<0 ||i1>NY-1 || j1 > NX-1 || k1>NZ-1){
      cDose3dr[tid]=0.0;
    } else {
      cx0=cay[i0];cx1=cay[i1];cy0=cax[j0];cy1=cax[j1];cz0=caz[k0];cz1=caz[k1];
      xc=(caay[tid]-cx0)/(cx1-cx0);
      yc=(caax[tid]-cy0)/(cy1-cy0);
      zc=(caaz[tid]-cz0)/(cz1-cz0);
      c00=cDose3d[i0*NX*NZ+j0*NZ+k0]*(1-xc)+cDose3d[i1*NX*NZ+j0*NZ+k0]*xc;
      c10=cDose3d[i0*NX*NZ+j1*NZ+k0]*(1-xc)+cDose3d[i1*NX*NZ+j1*NZ+k0]*xc;
      c01=cDose3d[i0*NX*NZ+j0*NZ+k1]*(1-xc)+cDose3d[i1*NX*NZ+j0*NZ+k1]*xc;
      c11=cDose3d[i0*NX*NZ+j1*NZ+k1]*(1-xc)+cDose3d[i1*NX*NZ+j1*NZ+k1]*xc;
      c0=c00*(1-yc)+c10*yc;
      c1=c01*(1-yc)+c11*yc;
      cDose3dr[tid]=c0*(1-zc)+c1*zc;
    }
    tid+=blockDim.x*gridDim.x;
  }
}
//
void interps(int npnts, float *caax, float *caay, float *caaz,float *cDose3dr, \
	     int NX, float *cax, int NY, float *cay,			\
	     int NZ, float *caz, float * cDose3d,float * mnmxs){
  float *d_caax, *d_caay, *d_caaz, *d_cdos3r;
  float *d_cax, *d_cay, *d_caz, *d_cdos3;
  float *d_mnmx;
  //float *h_cdx, *h_cdy, *h_cdz, *h_cdos2;
  //float *h_rax, *h_ray, *h_raz, *h_cdos3;
  //int i;
  cudaMalloc( (void**)&d_caax, npnts*sizeof(float));
  cudaMalloc( (void**)&d_caay, npnts*sizeof(float));
  cudaMalloc( (void**)&d_caaz, npnts*sizeof(float));
  cudaMalloc( (void**)&d_cdos3r,npnts*sizeof(float));
  cudaMalloc( (void**)&d_cax, NX*sizeof(float));
  cudaMalloc( (void**)&d_cay, NY*sizeof(float));
  cudaMalloc( (void**)&d_caz, NZ*sizeof(float));
  cudaMalloc( (void**)&d_cdos3,NZ*NX*NY*sizeof(float));
  cudaMalloc( (void**)&d_mnmx,6*sizeof(float));
  /*cudaHostAlloc( (void**)&h_cdx, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_cdy, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_cdz, npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_cdos2,npnts*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_rax, rx*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_ray, ry*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_raz, rz*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc( (void**)&h_cdos3,rz*rx*ry*sizeof(float),cudaHostAllocWriteCombined | cudaHostAllocMapped);
  for(i=0;i<npnts;i++){
    h_cdx[i]=cdx[i];
    h_cdy[i]=cdy[i];
    h_cdz[i]=cdz[i];
    h_cdos2[i]=cdos2d[i];
  }
  for(i=0;i<rx;i++)h_rax[i]=rax*/
  
  //copy to gpu
  cudaMemcpy(d_caax,caax,npnts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_caay,caay,npnts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_caaz,caaz,npnts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_cdos3r,cDose3dr,npnts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_cax,cax,NX*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_cay,cay,NY*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_caz,caz,NZ*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_cdos3,cDose3d,NX*NY*NZ*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_mnmx,mnmxs,6*sizeof(float),cudaMemcpyHostToDevice);
  
  

  interp<<<512,512>>>(npnts,d_caax,d_caay,d_caaz,d_cdos3r, \
		      NX,d_cax,NY,d_cay,NZ,d_caz,d_cdos3,d_mnmx);

  cudaMemcpy(cDose3dr,d_cdos3r,npnts*sizeof(float),cudaMemcpyDeviceToHost);
  //cudaDeviceReset();
  cudaFree(d_caax);cudaFree(d_caay);cudaFree(d_caaz);cudaFree(d_cdos3r);
  cudaFree(d_cax);cudaFree(d_cay);cudaFree(d_caz);cudaFree(d_cdos3);
  cudaFree(d_mnmx);
}

