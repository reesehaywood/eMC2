#include <stdio.h>
#include <sys/time.h>
/*get random integers for seeding the random number gen*/
unsigned long int random_seed()
{
  unsigned int seed;
  struct timeval tv;
  //FILE *devrandom;
  //if ((devrandom=fopen("/dev/random","r"))==NULL){
    gettimeofday(&tv,0);
    seed=tv.tv_sec+tv.tv_usec;
    //printf("Got seed %u from gettimeofday()\n",seed);
    //}else{
    //fread(&seed,sizeof(seed),1,devrandom);
    //printf("Got seed %u from /dev/random\n",seed);
    //fclose(devrandom);
    //}
  return(seed);
}
//point inside of a polygon for putting electrons in the block
__device__ int pnpoly(const int nvert, const float * __restrict__ vertx, const float * __restrict__ verty, const float testx, const float testy){
  int i,j,c=0;
  for(i=0,j=nvert-1;i<nvert;j=i++){
    if ( ( (verty[i]>testy) != (verty[j]>testy))&&(testx<(vertx[j]-vertx[i])*(testy-verty[i])/(verty[j]-verty[i])+vertx[i]) ){
      c= !c;
    }
  }
  return c;
}
//scale sesc to diff materials
__device__ float scalSeSc(const float e0, const int dens){
  /*float fdens=(float)dens;
int a=signbit(fdens-0.f)+signbit(0.0f-fdens);
int w=signbit(fdens-4.f)+signbit(4.0f-fdens);
int il=signbit(fdens-1.f)+signbit(1.f-fdens);
int el=signbit(fdens-2.f)+signbit(2.f-fdens);
int f=signbit(fdens-3.f)+signbit(3.f-fdens);
int m=signbit(fdens-5.f)+signbit(5.f-fdens);
int cb=signbit(fdens-6.f)+signbit(6.f-fdens);
int sb=signbit(fdens-7.f)+signbit(7.f-fdens);
int l=signbit(fdens-8.f)+signbit(8.f-fdens);
return (1.0f-(float)a)*(-0.0000008f*e0*e0*e0*e0+0.0000623f*e0*e0*e0-0.0017838f*e0*e0+0.0254018f*e0+0.8777285f)+(1.0f-(float)w)*1.0f+ \
       (1.0f-(float)il)*0.98907f+(1.0f-(float)el)*0.98f+ \
       (1.0f-(float)f)*(-0.00369f*__logf(e0)+1.02099f)+ \
       (1.0f-(float)m)*0.98872f+(1.0f-(float)cb)*(0.00463f*__logf(e0)+0.90144f)+ \
       (1.0f-(float)sb)*(0.00120f*__logf(e0)+0.92755f);
//not using lead*/
    switch(dens){
  case 0:
    //air
    return -0.0000008f*e0*e0*e0*e0+0.0000623f*e0*e0*e0-0.0017838f*e0*e0+0.0254018f*e0+0.8777285f;
  case 1:
    //inhaled lung
    return 0.98907f;
  case 2:
    //exhaled lung
    return 0.98907f;
  case 3:
    //fat
    return -0.00369f*__logf(e0)+1.02099f;
  case 4:
    //water
    return 1.0f;
  case 5:
    //muscle
    return 0.98872f;
  case 6:
    // cortical bone
    return 0.00463f*__logf(e0)+0.90144f;
  case 7:
    //skeletal bone
    return 0.00120f*__logf(e0)+0.92755f;
  case 8:
    //lead
    return 0.3056f*__logf(e0)+0.53970f;
  default:
    //return water
    return 1.0f;
    }
}
//scale sesy to diff materials
__device__ float scalSeSy(const float e0, const int dens){
  /*float fdens=(float)dens;
int a=signbit(fdens-0.f)+signbit(0.0f-fdens);
int w=signbit(fdens-4.f)+signbit(4.0f-fdens);
int il=signbit(fdens-1.f)+signbit(1.f-fdens);
int el=signbit(fdens-2.f)+signbit(2.f-fdens);
int f=signbit(fdens-3.f)+signbit(3.f-fdens);
int m=signbit(fdens-5.f)+signbit(5.f-fdens);
int cb=signbit(fdens-6.f)+signbit(6.f-fdens);
int sb=signbit(fdens-7.f)+signbit(7.f-fdens);
int l=signbit(fdens-8.f)+signbit(8.f-fdens);
return (1.0f-(float)a)*(-0.0000008f*e0*e0*e0*e0+0.0000623f*e0*e0*e0-0.0017838f*e0*e0+0.0254018f*e0+0.8777285f)+(1.0f-(float)w)*1.0f+ \
       (1.0f-(float)il)*1.011050785f+(1.0f-(float)el)*1.011050785f+ \
       (1.0f-(float)f)*(0.010376f*__logf(e0)+0.816777f)+ \
       (1.0f-(float)m)*1.01140869f+(1.0f-(float)cb)*(__fdividef(1.0f,(0.00463f*__logf(e0)+0.90144f)))+ \
       (1.0f-(float)sb)*(-0.0126f*__logf(e0)+1.342074f);
//not using lead*/
  switch(dens){
  case 0:
    //air
    return -0.0000008f*e0*e0*e0*e0+0.0000623f*e0*e0*e0-0.0017838f*e0*e0+0.0254018f*e0+0.8777285f;
  case 4:
    //water
    return 1.00f;
  case 1:
    //inhaled lung
    return 1.011050785f;
  case 2:
    //exhaled lung
    return 1.011050785f;
  case 3:
    //fat
     return 0.010376f*__logf(e0)+0.816777f;
  case 5:
    //muscle
    return 1.01140869f;
  case 6:
    //cortical bone
    return __fdividef(1.0f,(0.00463f*__logf(e0)+0.90144f));
  case 7:
    //skeletal bone
    return -0.0126f*__logf(e0)+1.342074f;
  case 8:
    //lead
    return -1.338112f*__logf(e0)+14.701651f;
  default:
    //return water if nothing found
    return 1.00f;
    }
}
///
//scale sesr to diff materials
__device__ float scalSeSr(const float e0, const int dens){
  /*float fdens=(float)dens;
int a=signbit(fdens-0.f)+signbit(0.0f-fdens);
int w=signbit(fdens-4.f)+signbit(4.0f-fdens);
int il=signbit(fdens-1.f)+signbit(1.f-fdens);
int el=signbit(fdens-2.f)+signbit(2.f-fdens);
int f=signbit(fdens-3.f)+signbit(3.f-fdens);
int m=signbit(fdens-5.f)+signbit(5.f-fdens);
int cb=signbit(fdens-6.f)+signbit(6.f-fdens);
int sb=signbit(fdens-7.f)+signbit(7.f-fdens);
int l=signbit(fdens-8.f)+signbit(8.f-fdens);
return (1.0f-(float)a)*(-0.0004f*e0+0.99583f)+(1.0f-(float)w)*1.0f+ \
       (1.0f-(float)il)*0.98792f+(1.0f-(float)el)*0.98792f+ \
       (1.0f-(float)f)*(0.00664f*__logf(e0)+0.83435f)+ \
       (1.0f-(float)m)*0.98620f+(1.0f-(float)cb)*(-0.01908f*__logf(e0)+1.42207f)+ \
       (1.0f-(float)sb)*(-0.01158f*__logf(e0)+1.24946f);*/
       //not using lead
    switch(dens){
  case 0:
    //air
    return -0.0004f*e0+0.99583f;
  case 4:
    //water
    return 1.00f;
  case 1:
    //inhaled lung
    return 0.98792f;
  case 2:
    //exhaled lung
    return 0.98792f;
  case 3:
    //fat
    return 0.00664f*__logf(e0)+0.83435f;
  case 5:
    //muscle
    return 0.98620f;
  case 6:
    //cortical bone
    return -0.01908f*__logf(e0)+1.42207f;
  case 7:
    //skeletal bone
    return -0.01158f*__logf(e0)+1.24946f;
  case 8:
    //lead
    if(e0<0.35f){
      return 17.85446f*e0+7.41971f;
    }else{
      return float(9.8963f*__powf(float(e0),-0.1661f));
    }
  default:
    //return water if nothing found
    return 1.00f;
    }
}

int wherb(const float *__restrict__ x,const float val, const int n){
  int nabove=n+1,nbelow=0,mid;
  while(nabove -nbelow >1){
    mid=(nabove+nbelow)/2;
    if(val==x[mid-1])return mid-1;
    if(val< x[mid-1]){ nabove=mid;}
    else{ nbelow=mid;}
  }
  return nbelow-1;
}
//utility to find an index in an array 1D on the GPU
//look up interpolations search in wikipedia it may be faster still
//binary search
__device__ int d_wherb(const float *__restrict__ x,const float val, const int n, int ind){
  int nabove=n+1,nbelow=0,mid;
  //pass in dlo/2
  //if(x[ind]<=val && val <x[ind+1]){
  //  return ind;
  //}else{
    while(nabove -nbelow >1){
      mid=(nabove+nbelow)/2;
      if(val==x[mid-1])return mid-1;
      if(val< x[mid-1]){ nabove=mid;}
      else{ nbelow=mid;}
    }
    return nbelow-1;
  //}
}
/*function defined by fitting the NIST table for radiative stopping power
  with a power log the fitting is done in Mathematica and has the form
  a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8+jx^9 
  where x is log(energy) or log(csda) as below
 */
__device__ float SeSr(const float sEn){
  float a=-4.354893286051908f,b=0.9401293148549328f;
  float c=0.1730498753429604f,d=-0.029628222973712436f;
  float e=-0.005884969822153593f,f=0.001492746368448323f;
  float g=-0.0001105126090672895f,h=-0.000034569280918223715f;
  float i=-2.3733670437520215e-7f,j=2.165906072712042e-7f;
  float lgx=__logf(sEn);
  float pgx2=lgx*lgx,pgx3=lgx*lgx*lgx;
  return __expf(a+__fmul_rn(b,lgx)+__fmul_rn(__fmul_rn(c,lgx),lgx)+ \
		__fmul_rn(d,pgx3)+__fmul_rn(e,pgx2*pgx2)+ \
		__fmul_rn(f,pgx3*pgx2)+__fmul_rn(g,pgx3*pgx3)+ \
		__fmul_rn(h,pgx2*pgx2*pgx3)+__fmul_rn(i,pgx3*pgx3*pgx2)+ \
		__fmul_rn(j,pgx3*pgx3*pgx3));
}

/*same as above for collisional stopping power*/
__device__ float SeSc(const float sEn){
  float a=0.6140175073693391f,b=-0.07427360869425056f;
  float c=0.08614598065631322f,d=-0.01995297887839599f;
  float e=-0.0006693170690464525f,f=0.0008310100951044747f;
  float g=-0.00005760622637751102f,h=-0.000012308838910962778f;
  float i=1.8037075818809332e-6f,j=-6.02581701567041e-8f;
  float lgx=__logf(sEn);
  float pgx2=lgx*lgx,pgx3=lgx*lgx*lgx;
  return __expf(a+__fmul_rn(b,lgx)+__fmul_rn(__fmul_rn(c,lgx),lgx)+ \
		__fmul_rn(d,pgx3)+__fmul_rn(e,pgx2*pgx2)+ \
		__fmul_rn(f,pgx3*pgx2)+__fmul_rn(g,pgx3*pgx3)+ \
		__fmul_rn(h,pgx2*pgx2*pgx3)+__fmul_rn(i,pgx3*pgx3*pgx2)+ \
		__fmul_rn(j,pgx3*pgx3*pgx3));
}
/*same as above for continuously slowing down approximation*/
__device__ float SeScsda(const float sEn){
  float a=-0.8286695741151627f,b=1.2334340420683476f;
  float c=-0.1038066835294922f,d=0.009672830002845856f;
  float e=0.0023275003264604252f,f=-0.0005957074416801652f;
  float g=-0.0000509683968743365f,h=0.000013328076877645397f;
  float i=3.2906109220493454e-7f,j=-9.558213646863889e-8f;
  float lgx=__logf(sEn);
  float pgx2=lgx*lgx,pgx3=lgx*lgx*lgx;
  return __expf(a+__fmul_rn(b,lgx)+__fmul_rn(__fmul_rn(c,lgx),lgx)+ \
		__fmul_rn(d,pgx3)+__fmul_rn(e,pgx2*pgx2)+ \
		__fmul_rn(f,pgx3*pgx2)+__fmul_rn(g,pgx3*pgx3)+ \
		__fmul_rn(h,pgx2*pgx2*pgx3)+__fmul_rn(i,pgx3*pgx3*pgx2)+ \
		__fmul_rn(j,pgx3*pgx3*pgx3));
}
/*same as above for radiation percent yield*/
__device__ float SeSy(const float sEn){
  float a=-5.627326516537827f,b=0.9207796926434323f;
  float c=0.08451335637434465f,d=0.000686569407791968f;
  float e=-0.006080923668522947f,f=-0.00006121087489969487f;
  float g=0.00021291343326992302f,h=-8.629294419917445e-6f;
  float i=-3.3109708642557386e-6f,j=2.813742105966361e-7f;
  float lgx=__logf(sEn);
  float pgx2=lgx*lgx,pgx3=lgx*lgx*lgx;
  return __expf(a+__fmul_rn(b,lgx)+__fmul_rn(__fmul_rn(c,lgx),lgx)+ \
		__fmul_rn(d,pgx3)+__fmul_rn(e,pgx2*pgx2)+ \
		__fmul_rn(f,pgx3*pgx2)+__fmul_rn(g,pgx3*pgx3)+ \
		__fmul_rn(h,pgx2*pgx2*pgx3)+__fmul_rn(i,pgx3*pgx3*pgx2)+ \
		__fmul_rn(j,pgx3*pgx3*pgx3));
}
/*fit of energy vs CSDA i.e. the inverse of the one above*/
__device__ float ScsdaSe(const float sCSDA){
  float a=0.7123832103497709f,b=0.9052245441729793f;
  float c=0.05649056449346559f,d=-0.000644522927760566f;
  float e=-0.0009342635644128602f,f=0.0000679237444644867f;
  float g=0.00005598568823350923f,h=9.826586704508063e-6f;
  float i=8.371413097783667e-7f,j=2.9186499738516662e-8f;
  float lgx=__logf(sCSDA);
  float pgx2=lgx*lgx,pgx3=lgx*lgx*lgx;
  return __expf(a+__fmul_rn(b,lgx)+__fmul_rn(__fmul_rn(c,lgx),lgx)+ \
		__fmul_rn(d,pgx3)+__fmul_rn(e,pgx2*pgx2)+ \
		__fmul_rn(f,pgx3*pgx2)+__fmul_rn(g,pgx3*pgx3)+ \
		__fmul_rn(h,pgx2*pgx2*pgx3)+__fmul_rn(i,pgx3*pgx3*pgx2)+ \
		__fmul_rn(j,pgx3*pgx3*pgx3));
}

/*utility to find the max value in an array 1D*/
float maxval(const float x[],const int n){
  float max;
  int i;
  max=x[0];
  for (i=0;i<n;i++){
    if(x[i]>=max){
      max=x[i];
    }
  }
  return(max);
}
