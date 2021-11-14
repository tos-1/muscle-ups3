#include <stdio.h>
#include <stdlib.h>
#include "pnum.h"

#ifndef pbc
#define pbc(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))
#endif

/* https://www.tutorialspoint.com/computer_graphics/circle_generation_algorithm.htm */
int tagDisk( int NG, int xyz_c[], int r, int *sph, int id, int targid ){
  int x = 0, y = r;
  int d = 3 - 2 * r;
  int i, index;
  int collapse = targid;
  int M=0;	

  i=-abs(x);
  while(i<=abs(x)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+y,NG),xyz_c[2],NG);
    if (*(sph+index)==collapse){
      *(sph+index)=id;
      M++;
    }
    i++;
  }
  i=-abs(x);
  while(i<=abs(x)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-y,NG),xyz_c[2],NG);
    if (*(sph+index)==collapse){
      *(sph+index)=id;
      M++;
    }
    i++;
  }
  i=-abs(y);
  while(i<=abs(y)){
    index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+x,NG),xyz_c[2],NG);
    if (*(sph+index)==collapse){
      *(sph+index)=id;
      M++;
    }
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
      if (*(sph+index)==collapse){
        *(sph+index)=id;
        M++;
      }
      i++;
    }
    i=-abs(x);
    while(i<=abs(x)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-y,NG),xyz_c[2],NG);
      if (*(sph+index)==collapse){
        *(sph+index)=id;
        M++;
      }
      i++;
    }
    i=-abs(y);
    while(i<=abs(y)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]+x,NG),xyz_c[2],NG);
      if (*(sph+index)==collapse){
        *(sph+index)=id;
        M++;
      }
      i++;
    }
    i=-abs(y);
    while(i<=abs(y)){
      index = pnum(pbc(xyz_c[0]+i,NG),pbc(xyz_c[1]-x,NG),xyz_c[2],NG);
      if (*(sph+index)==collapse){
        *(sph+index)=id;
        M++;
      }
      i++;
    }
  }
  return M;
}


/* Function for semi-circle generation, using Bresenham's algorithm */
int tagSphere( int NG, int index, int r , int *sph, int id, int targid ){
  int M=0;
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
  M+=tagDisk(NG, xyz_c2, x, sph, id, targid);
  xyz_c2[2] = pbc(xyz_c[2]-z,NG);
  M+=tagDisk(NG, xyz_c2, x, sph, id, targid);
  xyz_c2[2] = pbc(xyz_c[2]+x,NG);
  M+=tagDisk(NG, xyz_c2, z, sph, id, targid);

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
    M+=tagDisk(NG, xyz_c2, x, sph, id, targid);
    xyz_c2[2] = pbc(xyz_c[2]-z,NG);
    M+=tagDisk(NG, xyz_c2, x, sph, id, targid);
    xyz_c2[2] = pbc(xyz_c[2]+x,NG);
    M+=tagDisk(NG, xyz_c2, z, sph, id, targid);
    xyz_c2[2] = pbc(xyz_c[2]-x,NG);
    M+=tagDisk(NG, xyz_c2, z, sph, id, targid);
  }
  return M;
}
