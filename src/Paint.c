#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pbc(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))

int pnum(int i, int j, int k, int NG){
  return NG*NG*i + NG*j + k;
}


void Disk( int NG, int xyz_c[], int r, float *sph){
  int x = 0, y = r;
  int d = 3 - 2 * r;
  int i, index;
  float collapse = -3.0;
  
  i=-abs(x);
  while(i<=abs(x)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+y,NG),xyz_c[2],NG);
    *(sph+index)=collapse;
    i++;
  }
  i=-abs(x);
  while(i<=abs(x)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-y,NG),xyz_c[2],NG);
    *(sph+index)=collapse;
    i++;
  }
  i=-abs(y);
  while(i<=abs(y)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+x,NG),xyz_c[2],NG);
    *(sph+index)=collapse;
    i++;
  }

  while (y >= x){
    x++;
  
    if (d > 0) {
      y--;
      d = d + 4 * (x - y) + 10;
    }
    else {
      d = d + 4 * x + 6;
    }
    i=-abs(x);
    while(i<=abs(x)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+y,NG),xyz_c[2],NG);
      *(sph+index)=collapse;
      i++;
    }
    i=-abs(x);
    while(i<=abs(x)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-y,NG),xyz_c[2],NG);
      *(sph+index)=collapse;
      i++;
    }
    i=-abs(y);
    while(i<=abs(y)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+x,NG),xyz_c[2],NG);
      *(sph+index)=collapse;
      i++;
    }
    i=-abs(y);
    while(i<=abs(y)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-x,NG),xyz_c[2],NG);
    *(sph+index)=collapse;
    i++;
    }
  }
}


/* Function for semi-circle generation, using Bresenham's algorithm */ 
void Sphere( int NG, int index, int r , float *sph ){
  int xyz_c[3],xyz_c2[3];
  xyz_c[2] = index%NG;
  xyz_c[1] = ((index-xyz_c[2])/NG)%NG;
  xyz_c[0] = ( index-xyz_c[1]*NG-xyz_c[2] )/NG/NG;

  xyz_c2[0] = xyz_c[0];
  xyz_c2[1] = xyz_c[1];

  /* fixed yc */
  int x = 0, z = r;
  int d = 3 - 2 * r;

  /* don't want to miss the starting points */
  xyz_c2[2] = pbc(xyz_c[2]+z,NG);
  Disk(NG, xyz_c2, x, sph);
  xyz_c2[2] = pbc(xyz_c[2]-z,NG);
  Disk(NG, xyz_c2, x, sph);
  xyz_c2[2] = pbc(xyz_c[2]+x,NG);
  Disk(NG, xyz_c2, z, sph);

  while (z >= x){
    x++; 
    
    if (d > 0){
      z--;
      d = d + 4 * (x - z) + 10;
    }
    else{
      d = d + 4 * x + 6;
    }
    xyz_c2[2] = pbc(xyz_c[2]+z,NG);
    Disk(NG, xyz_c2, x, sph);
    xyz_c2[2] = pbc(xyz_c[2]-z,NG);
    Disk(NG, xyz_c2, x, sph);
    xyz_c2[2] = pbc(xyz_c[2]+x,NG);
    Disk(NG, xyz_c2, z, sph);
    xyz_c2[2] = pbc(xyz_c[2]-x,NG);
    Disk(NG, xyz_c2, z, sph);
  }
}


int Paint( int NG, int ns, int *sift, float *psi ){
  int NGt3 = NG*NG*NG;
  int i;
  int R;

  for (i=0;i<NGt3;i++){
    R = *(sift+i);
    if (R>=ns){
      Sphere( NG, i, R, psi );
    }
  }
  printf("finished painting halos\n");
  return 0;
}
