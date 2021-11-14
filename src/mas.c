#include <math.h>
#include <stdio.h>
#include "mas.h"
#define min(A,B) (((A)<(B)) ? (A):(B))

void mas( int ng, int np, double boxsize, float *pos, float *density){
  float cell_len= (float) boxsize/ng;
  float dngx,dngy,dngz;
  int NP = np*np*np;   //number of particles
  int NG = ng*ng*ng;   //number of grid cells
  float ratio = (float)NG/NP;
  float x,y,z;
  float minden=0.;
  int cell[3];
  int index, i, count=0;

  for (int i=0; i<NP; i++){

    x = *(pos + i );
    y = *(pos + NP + i );
    z = *(pos + 2*NP + i );

    cell[0]=(int) floor(x/cell_len);
    cell[1]=(int) floor(y/cell_len);
    cell[2]=(int) floor(z/cell_len);

    dngx = x/(cell_len) - (float)cell[0];
    dngy = y/(cell_len) - (float)cell[1];
    dngz = z/(cell_len) - (float)cell[2];

    index = cell[2] + ng * cell[1] + ng*ng* cell[0];
    *(density+index) += ratio*(1-dngx)*(1-dngy)*(1-dngz);
    index = cell[2] + ng * cell[1] + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += ratio*dngx*(1-dngy)*(1-dngz);
    index = cell[2] + ng *((cell[1]+1)%ng) + ng*ng* cell[0];
    *(density+index) += ratio*(1-dngx)*dngy*(1-dngz);
    index = ((cell[2]+1)%ng) + ng * cell[1] + ng*ng* cell[0];
    *(density+index) += ratio*(1-dngx)*(1-dngy)*dngz;
    index = cell[2] + ng * ((cell[1]+1)%ng) + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += ratio*dngx*dngy*(1-dngz);
    index = ((cell[2]+1)%ng) + ng * cell[1] + ng*ng* ((cell[0]+1)%ng);
    *(density + index) += ratio*dngx*(1-dngy)*dngz;
    index = ((cell[2]+1)%ng) + ng * ((cell[1]+1)%ng) + ng*ng* cell[0];
    *(density+ index) += ratio*(1-dngx)*dngy*dngz;
    index = ((cell[2]+1)%ng) + ng * ((cell[1]+1)%ng) + ng*ng* ((cell[0]+1)%ng);
    *(density+index) += ratio*dngx*dngy*dngz;
  }
  // get overdensity
  for (i=0; i<NG; i++){
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
    if (*(density+i)<=-1.0) *(density+i)=minden;
  }
  printf("ratio of particles dens<-1 is %f\n",(float)count/NP);
}
