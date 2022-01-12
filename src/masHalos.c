#include <math.h>
#include <stdio.h>
#include "masHalos.h"
#define min(A,B) (((A)<(B)) ? (A):(B))

void masHalos( int ng, int nH, double boxsize, float *pos, float *density, int *Nx){
  float cell_len= (float) boxsize/(float) ng;
  float dngx,dngy,dngz;
  int NG = ng*ng*ng;   //number of grid cells
  float x,y,z;
  float minden=0.;
  int cell[3];
  int index, i, count=0;
  int nx, avgnx=0;
  float favgnx;

  for (int i=0; i<nH; i++){

    nx = *(Nx+i);
    avgnx += nx;

    x = *(pos + i );
    y = *(pos + nH + i );
    z = *(pos + 2*nH + i );

    cell[0]=(int) floor(x/cell_len);
    cell[1]=(int) floor(y/cell_len);
    cell[2]=(int) floor(z/cell_len);

    dngx = x/(cell_len) - (float)cell[0];
    dngy = y/(cell_len) - (float)cell[1];
    dngz = z/(cell_len) - (float)cell[2];

    index = cell[2] + ng * cell[1] + ng*ng* cell[0];
    *(density+index) += nx * (1-dngx)*(1-dngy)*(1-dngz);
    index = cell[2] + ng * cell[1] + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += nx * dngx*(1-dngy)*(1-dngz);
    index = cell[2] + ng *((cell[1]+1)%ng) + ng*ng* cell[0];
    *(density+index) += nx * (1-dngx)*dngy*(1-dngz);
    index = ((cell[2]+1)%ng) + ng * cell[1] + ng*ng* cell[0];
    *(density+index) += nx * (1-dngx)*(1-dngy)*dngz;
    index = cell[2] + ng * ((cell[1]+1)%ng) + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += nx * dngx*dngy*(1-dngz);
    index = ((cell[2]+1)%ng) + ng * cell[1] + ng*ng* ((cell[0]+1)%ng);
    *(density + index) += nx * dngx*(1-dngy)*dngz;
    index = ((cell[2]+1)%ng) + ng * ((cell[1]+1)%ng) + ng*ng* cell[0];
    *(density+ index) += nx * (1-dngx)*dngy*dngz;
    index = ((cell[2]+1)%ng) + ng * ((cell[1]+1)%ng) + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += nx * dngx*dngy*dngz;
  }

  // get average number of halo particles
  favgnx = (float) avgnx/nH;
  printf("average number of halo particles %f\n",favgnx);

  // normalize
  for (i=0; i<NG; i++){
    *(density+i)*=(float)NG/favgnx/nH;
    *(density+i)-=1.0;
    if (*(density+i)>-1.0){
      minden=min(*(density+i),minden);	
      continue;
    }
    else{
      count++;
    }
  }
  printf("minden=%f\n",minden);
  for (i=0; i<NG; i++){
    //printf("%f\n",*(density+i));
    if (*(density+i)<-1.0) *(density+i)=minden;
  }
  printf("ratio of particles with density<-1 is %f\n",(float) count/nH);
}
