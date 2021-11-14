#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include "pnum.h"
#include "search.h"
#include "sort.h"
#include "bresenham.h"

#define pbc(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))
#define diffpbc(A,B) ( ( (fabs(A)) < (((B)/(2))) ) ? (A):( (A>0) ? (fabs(A)-B):(B-fabs(A)) ) )
#define SQRT(A) sqrt(A)
#define SIN(A) sin(A)
#define COS(A) cos(A)
#define ACOS(A) acos(A)
#define ATAN(A,B) atan2(A,B)
#define PI 3.14159265358979323846
#define MHP 20


typedef struct{
    int size;
    int id;
    int npinh;
    int count;
    double Mass;
    double Rv;
    double cM;
    double poshc[3];
    int halopid[];
} Halo;

typedef struct{
    int np;
    int count;
    int *pid;
    float *cc;
} SameSmooth;


int halomodel( const int NG, const double boxsize, const double redshift, const double Delta_v, const double mcrit, const double rho, 
		int *sift, int *psi3, float *cc, double *pos, const char *path, const char *hmf){
  int nhp, tnh, tnhs, nres, npinh;
  int i, j, k, h, l, m, idcan;
  int *halonum, *numinh;
  int *tmp_countsm, *numinh_tmp;
  int *haloID, *halop;
  double poshcx , poshcy , poshcz;
  double Dx, Dy, Dz;
  double phi, theta;
  double q, r, Rvir, rk, M1;
  double x , mass, cm , y;
  double *dist0c;
  double dnpinh;
  double xyz[3];
  //double axisx,axisy,axisz;
  //double Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Iij[9];
  //double evals[3], mevac[9], xyz2[3];
  int *id0c,status;
  char saveToPath_halocatalogue[500]={'\0'};
  char saveToPath_halonum[500]={'\0'};
  const int nothalonum = NG*NG*NG;
  double dnothalonum = (double) nothalonum;
  const double unit_mass = rho * boxsize*boxsize*boxsize/dnothalonum;
  printf("unit_mass=%.4e\n",unit_mass);

  /* count halo particles per smoothing */
  nhp = 0;
  tmp_countsm = malloc(100*sizeof(int *));

  if (tmp_countsm==NULL){
    printf("Failed to allocate memory for tmp_countsm.");
    exit(0);
  }

  for (i=0;i<100;i++){
    tmp_countsm[i]=0;
  }
  for (i=0;i<nothalonum;i++){
    if (*(sift + i) != -1){
      nhp ++;
      tmp_countsm[*(sift + i)] += 1;
    }
  }
  printf("number of halo seeds=%d\n", nhp);
  tnhs=0; //# of smoothing scales
  for (i=0;i<100;i++){
    if (tmp_countsm[i]!=0){
      tnhs++;
    }
  }
  const int nhs = tnhs;
  printf("number of smoothings found=%d\n", nhs);
  int countsm[nhs];//# of candidate halos per smoothing scale
  j=0;
  for (i=0;i<100;i++){
    if (tmp_countsm[i]!=0){
      countsm[j]=tmp_countsm[i];
      tmp_countsm[i] = j;//to use later as bookmark for psm
      j++;
    }
  }

  /* store halo seeds in array */
  halop = malloc(nhp*sizeof(int *));
  if (halop==NULL){
    printf("Failed to allocate memory for halop.");
    exit(0);
  }

  j = 0;
  for (i=0;i<nothalonum;i++){
    if (*(sift + i) != -1){
      *(halop+j) = i;
      j++;
    }
  }

  if (nhp!=j){
    printf("nhp and j must be equal");
    exit(0);
  }

  /* store halo seeds based on smoothing */
  //array of nhs pointers to struct SameSmooth
  SameSmooth *sm[nhs]; 
  for (h=0;h<nhs;h++){
    sm[h] = malloc( sizeof(SameSmooth) );
    if (sm[h]==NULL){
      printf("Failed to allocate memory for sm[h].");
      exit(0);
    }
  }

  //pointer to array of pointers to struct sm 
  SameSmooth *(*psm)[] = &sm;
  for (h=0;h<nhs;h++){
    k = countsm[h];
    (*psm)[h]->np=k;
    (*psm)[h]->count = 0;
    (*psm)[h]->pid = malloc( k*sizeof(int *) );
    (*psm)[h]->cc = malloc( k*sizeof(float *) );
    if ((*psm)[h]->pid==NULL){
      printf("Failed to allocate memory for (*psm)[h]->pid.");
      exit(0);
    }
    if ((*psm)[h]->cc==NULL){
      printf("Failed to allocate memory for (*psm)[h]->cc.");
      exit(0);
    }
    for (i=0;i<k;i++){
      (*psm)[h]->pid[i] = 0;
      (*psm)[h]->cc[i] = 0.0;
    }
  }

  for (h=0;h<nhp;h++){
    j = *(halop+h);
    k = tmp_countsm[*(sift+j)];
    i = (*psm)[k]->count;
    (*psm)[k]->pid[i]=j;
    (*psm)[k]->cc[i]=*(cc+j);
    (*psm)[k]->count++;
  }
  free(tmp_countsm);
  free(halop);

  for (h=0;h<nhs;h++){
    if ((*psm)[h]->count != (*psm)[h]->np){
      printf("You missed some particles\n");
      exit(0);
    }
  }

  /* Reorder them based on collapse criterion */
  for (h=0;h<nhs;h++){
    j = (*psm)[h]->np;
    mergeSortFloat( (*psm)[h]->cc, (*psm)[h]->pid, 0, j-1 );
  }

  /* read lookup halo mass function */
  int nlogMbins, *noflogM, *nM;
  double dlogM, *logM, *dn0dlnM;
  FILE *fp;
  fp = fopen(hmf,"rb");

  if (fp == NULL){
    printf("Could not open file\n");
    exit(0);
  }
  fread(&nlogMbins, sizeof(int), 1, fp);//this are actually bin edges
  nlogMbins-=1;
  if (nlogMbins<2){
    printf("select a reasonable number of bins of the HMF.\n");
  }

  logM = malloc((nlogMbins+1)*sizeof(double *));//bin edges
  if (logM==NULL){
    printf("Failed to allocate memory for logM.\n");
    exit(42);
  }

  dn0dlnM = malloc((nlogMbins+1)*sizeof(double *));//actual hmf
  if (dn0dlnM==NULL){
    printf("Failed to allocate memory dn0dlnM.\n");
    exit(42);
  }

  nM = malloc(nlogMbins*sizeof(int *)); //hmf in terms of particles
  if (nM==NULL){
    printf("Failed to allocate memory for nM.\n");
    exit(42);
  }

  noflogM = malloc((nlogMbins+1)*sizeof(int *));//bin edges in terms of particles
  if (noflogM==NULL){
    printf("Failed to allocate memory noflogM.\n");
    exit(42);
  }

  for (i=0;i<nlogMbins+1;i++){
    fread(logM+i, sizeof(double), 1, fp);
  }

  for (i=0;i<nlogMbins+1;i++){
    *(noflogM+i) = round(pow(10,*(logM+i))/unit_mass);
  }
  dlogM = *(logM+1)-*logM;
  double dlnM = log(pow(10,dlogM)); //natural delta log

  for (i=0;i<nlogMbins+1;i++){
    fread(dn0dlnM+i, sizeof(double), 1, fp);
  }
  fclose(fp);


  double tmp;
  for (i=0;i<nlogMbins;i++){
    tmp=(log10(dn0dlnM[i+1])-log10(dn0dlnM[i]))/(logM[i+1]-logM[i])*dlogM/2.+log10(dn0dlnM[i]);
    tmp=pow(10,tmp);
    *(nM+i)=round(tmp*dlnM*boxsize*boxsize*boxsize);
  }
  free(logM);
  free(dn0dlnM);

  /* halo-id array */
  j = 0;
  halonum = malloc(nothalonum*sizeof(int *));
  for (i=0;i<nothalonum;i++){
    if (*(psi3+i)==0){//no collapse
      *(halonum+i) = nothalonum;
    }
    else{//collapse
      *(halonum+i) = -1;
    }
  }

  int nhc=0;
  for (h=0;h<nhs;h++){
    nhc += (*psm)[h]->np;//tot # of halo candidates
  }

  int *mhc,*idhc,*bhnfc,*hnf;
  mhc = malloc( nhc*sizeof(int *) );//mass halo candidates
  if (mhc==NULL){
    printf("Failed to allocate memory for mhc.");
    exit(42);
  }

  idhc = malloc( nhc*sizeof(int *) );//id halo candidates
  if (idhc==NULL){
    printf("Failed to allocate memory for idhc.");
    exit(42);
  }

  bhnfc = malloc(nhc*sizeof(int *));
  if (bhnfc==NULL){
    printf("Failed to allocate memory for bhnfc.");
    exit(42);
  }

  hnf = malloc( nlogMbins*sizeof(int *) );//halo number function
  if (hnf==NULL){
    printf("Failed to allocate memory hnf.");
    exit(42);
  }

  for (h=0;h<nhc;h++){
    *(mhc+h)=0;
    *(idhc+h)=-1;
    *(bhnfc+h)=-1;
  }

  for (h=0;h<nlogMbins;h++){
    *(hnf+h)=0;
  }

  int *eps, eps1, eps0;
  eps = malloc(nlogMbins * sizeof(int *));

  if (eps==NULL){
    printf("Failed to allocate memory.");
    exit(42);
  }

  for (i=0;i<nlogMbins;i++){
      eps[i] = nM[i];
  }
  
  /* compile a catalogue of halos without mergings */
  /*tnh=0;
  int M;
  for (h=nhs-1;h>=0;h--){
    for ( i=0; i<(*psm)[h]->np; i++ ){//candidates
      idcan = (*psm)[h]->pid[i];
      if (*(halonum+idcan)==-1){
        M = tagSphere( NG, idcan, *(sift+idcan), halonum, idcan, -1);
	*(mhc + tnh )=M;
	*(idhc + tnh )=idcan;
        tnh++;
      }
      else{//merge if overlap
         M = tagSphere( NG, idcan, *(sift+idcan), halonum,  *(halonum+idcan), -1);
      }
    }
  }
  free(eps);
  free(bhnfc);
  free(hnf);
  free(noflogM);
  free(nM);
  const int nh = tnh;
  printf("halo catalogue compiled, number of halos=%d\n", nh);
  free(mhc);
  free(idhc);
  for (h=0;h<nhs;h++){
    free((*psm)[h]->pid);
    free((*psm)[h]->cc);
    free(sm[h]);
  }*/

  /* compile a catalogue */
  tnh=0;
  int M, ctnh, rs, nmerg=0;
  for (h=nhs-1;h>=0;h--){
    k=h;
    while(k>=0){
      l = (*psm)[k]->pid[0];
      rs = *( sift + l );
      for ( i=0; i<(*psm)[h]->np; i++){//candidates
        idcan = (*psm)[h]->pid[i];
 
        if (*(halonum+idcan)==-1){//host halo
          M = tagSphere( NG, idcan, *(sift+idcan), halonum, idcan, -1);
          *(mhc + tnh)=M;
          for (j=0;j<nlogMbins;j++){//update hnf
            if ( (M>=*(noflogM+j)) && (M<*(noflogM+j+1)) ){
              *(hnf+j)+=1;
              eps0=abs(*(nM+j)-*(hnf+j));
              break;
            }
            else if (M>=*(noflogM+nlogMbins)){//exceeds max bin
              j=nlogMbins-1;
              *(hnf+j)+=1;
              eps0=abs(*(nM+j)-*(hnf+j));
            }
          }
          if (eps0<=*(eps+j)){//add to catalogue
            *(idhc +tnh)=idcan;
            *(bhnfc+tnh)=j;
            tnh++;
            *(eps+j)=eps0;
            insertSort(idhc,mhc,bhnfc,tnh);//unlike mergeSort here goes n
          }
          else{//remove from catalogue
            M = tagSphere( NG, idcan, *(sift+idcan), halonum, -1, idcan);//untag
            *(mhc + tnh)=0;
            *(hnf+j)-=1;
          }
        }
  
        else if (*(sift+idcan)==rs){//merge?
          M = tagSphere( NG, idcan, *(sift+idcan), halonum, idcan, -1);//at first keep idcan
          ctnh = binarySearch( idhc, 0, tnh-1, *(halonum+idcan));
          *(mhc + ctnh)+=M;
          j=*(bhnfc+ctnh);
          for (m=0;m<nlogMbins;m++){//update hnf
            if ( (*(mhc + ctnh)>=*(noflogM+m)) && (*(mhc + ctnh)<*(noflogM+m+1)) ){
              break;
            }
          }
          if (j!=m){//changed bin
            *(hnf+j)-=1;
            *(hnf+m)+=1;
            eps0=abs(*(nM+j)-*(hnf+j));
            eps1=abs(*(nM+m)-*(hnf+m));
            if (eps1<=*(eps+m)){//
              M = tagSphere( NG, idcan, *(sift+idcan), halonum, *(halonum+idcan), idcan);//now tag with host id
              *(bhnfc+ctnh)=m;
              *(eps+j)=eps0;
              *(eps+m)=eps1;
              nmerg++;
            }
            else{
              M = tagSphere( NG, idcan, *(sift+idcan), halonum, -1, idcan);//untag
              *(mhc + ctnh)-=M;
              *(hnf+j)+=1;
              *(hnf+m)-=1;
            }
          }
          else{//did not change bin
            M = tagSphere( NG, idcan, *(sift+idcan), halonum, *(halonum+idcan), idcan);//now tag with host id
            nmerg++;
          }
        }
      }
      k--;
    }
  }
  free(eps);
  free(bhnfc);
  free(hnf);
  free(noflogM);
  free(nM);
  printf("number of mergings=%d\n", nmerg);
  const int nh = tnh;
  printf("halo catalogue compiled, number of halos=%d\n", nh);
  free(mhc);
  free(idhc);
  for (h=0;h<nhs;h++){
    free((*psm)[h]->pid);
    free((*psm)[h]->cc);
    free(sm[h]);
  }
  /* save halos ID and count halo particles */
  haloID = malloc( nh *sizeof(int *) );
  if (haloID==NULL){
    printf("Failed to allocate memory for haloID.");
    exit(0);
  }

  numinh = malloc( nh *sizeof(int *) );
  if (numinh==NULL){
    printf("Failed to allocate memory for numinh.");
    exit(0);
  }

  numinh_tmp = malloc( nothalonum *sizeof(int *) );
  if (numinh_tmp==NULL){
    printf("Failed to allocate memory for numinh_tmp.");
    exit(0);
  }

  for (i=0;i<nh;i++){
    haloID[i] = 0;
    numinh[i] = 0;
  }

  for (i=0;i<nothalonum;i++){
    numinh_tmp[i] = 0;
  }

  h = 0;
  nres = 0;
  for (i=0;i<nothalonum;i++){
    if ( (*(halonum+i) != nothalonum) && (*(halonum+i)!=-1) ){
      *(numinh_tmp+*(halonum+i)) = *(numinh_tmp+*(halonum+i)) + 1;
    }
    if (*(halonum+i) == i){
      *(haloID+h) = i;
      h++;
    }
    if (*(halonum+i) == -1){
      nres++;
    }
  }

  //if (nres!=0){
    //printf("halo particles not assigned to halos=%d\n", nres);
    //exit(0);
  //}

  if (h!=nh){
    printf("h and nh must be equal\n");
    exit(0);
  }

  for (h=0;h<nh;h++){
    *(numinh+h) = *(numinh_tmp + *(haloID+h));
  }
  free(numinh_tmp);

  /* sort haloID and numinh */
  mergeSortInt( haloID, numinh, 0, nh-1 );

  /* create a structure for halos */
  Halo *halos[nh];

  for (h=0;h<nh;h++){
    k = *(numinh+h);
    halos[h] = malloc( sizeof(Halo) + k*sizeof(int) );
    if (halos[h]==NULL){
      printf("Failed to allocate memory for halos[h].");
      exit(0);
    }
  }

  Halo *(*phalos)[] = &halos;
  for (h=0;h<nh;h++){
    k = *(numinh+h);
    (*phalos)[h]->size = sizeof(Halo) + k*sizeof(int);
    (*phalos)[h]->id=*(haloID+h);
    (*phalos)[h]->count=0;
    (*phalos)[h]->npinh=k;
    (*phalos)[h]->Mass=0.;
    (*phalos)[h]->Rv=0.;
    (*phalos)[h]->cM=0.;
    for (j=0;j<3;j++){
      (*phalos)[h]->poshc[j]=0.;
    }
    for (i=0;i<k;i++){
      (*phalos)[h]->halopid[i]=0;
    }
  }

  /* assign each halo particle to a halo */
  for (i=0;i<nothalonum;i++){
    if (*(halonum+i)!=nothalonum && (*(halonum+i)!=-1) ){
      h = binarySearch( haloID, 0, nh-1, *(halonum+i));
      (*phalos)[h]->halopid[(*phalos)[h]->count] = i;
      (*phalos)[h]->count++;
    }
  }

  /* check you counted them correctly */
  for (h=0;h<nh;h++){
    if ((*phalos)[h]->count != (*phalos)[h]->npinh){
      printf("You missed to add some particles in the halo\n");
      printf("Stopped at h=%d\n",h);
      exit(0);
    }
  }

  /* find the center of mass of each halo */
  for (h=0;h<nh;h++){
    npinh = (*phalos)[h]->npinh;
    dnpinh = (double) npinh;
    idcan = (*phalos)[h]->id;
    poshcx = *(pos+idcan);
    poshcy = *(pos+idcan+nothalonum);
    poshcz = *(pos+idcan+2*nothalonum);
    Dx=0., Dy=0., Dz=0.;
    for (i=0;i<npinh;i++){
      j = (*phalos)[h]->halopid[i];
      Dx+=diffpbc(*(pos+j)-poshcx,boxsize);
      Dy+=diffpbc(*(pos+j+nothalonum)-poshcy,boxsize);
      Dz+=diffpbc(*(pos+j+2*nothalonum)-poshcz,boxsize);
    }
    (*phalos)[h]->poshc[0]=pbc(poshcx+Dx/dnpinh,boxsize);
    (*phalos)[h]->poshc[1]=pbc(poshcy+Dy/dnpinh,boxsize);
    (*phalos)[h]->poshc[2]=pbc(poshcz+Dz/dnpinh,boxsize);
  }

  /* reorder particles according to NFW */
  for (h=0;h<nh;h++){
    npinh = (*phalos)[h]->npinh;
    dist0c = malloc(npinh*sizeof(double*));
    if (dist0c==NULL){
      printf("failed to allocate memory for dist0c\n");
      exit(42);
    }
    id0c = malloc(npinh*sizeof(int*));
    if (id0c==NULL){
      printf("failed to allocate memory for id0c\n");
      exit(42);
    }
    poshcx=(*phalos)[h]->poshc[0];
    poshcy=(*phalos)[h]->poshc[1];
    poshcz=(*phalos)[h]->poshc[2];

    //https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
    //Ixx=0.0, Iyy=0.0, Izz=0.0, Ixy=0.0, Ixz=0.0, Iyz=0.0;
    for (i=0;i<npinh;i++){
      j = (*phalos)[h]->halopid[i];
      Dx=*(pos+j)-poshcx;
      Dy=*(pos+j+nothalonum)-poshcy;
      Dz=*(pos+j+2*nothalonum)-poshcz;
      Dx=diffpbc(Dx,boxsize);
      Dy=diffpbc(Dy,boxsize);
      Dz=diffpbc(Dz,boxsize);
      *(id0c+i) = j;
      *(dist0c+i) = SQRT(Dx*Dx + Dy*Dy + Dz*Dz);
      //Ixx+=Dy*Dy+Dz*Dz, Iyy+=Dx*Dx+Dz*Dz, Izz+=Dx*Dx+Dy*Dy;
      //Ixy-=Dx*Dy, Ixz-=Dx*Dz, Iyz-=Dy*Dz;
    }
    mass = npinh * unit_mass;
    cm = 9./(1.+redshift)*pow(mass/mcrit, -0.13);
    cm*=2.0;
    Rvir = pow( 3. * mass / 4. / PI / Delta_v / rho, 1./3.);
    (*phalos)[h]->Mass=mass;
    (*phalos)[h]->Rv=Rvir;
    (*phalos)[h]->cM=cm;

    if (npinh<MHP) continue;//
    //Iij[0]=unit_mass*Ixx,Iij[1]=unit_mass*Ixy,Iij[2]=unit_mass*Ixz,Iij[3]=unit_mass*Ixy;Iij[4]=unit_mass*Iyy;
    //Iij[5]=unit_mass*Iyz,Iij[6]=unit_mass*Ixz,Iij[7]=unit_mass*Iyz,Iij[8]=unit_mass*Izz;

    //gsl_matrix_view I = gsl_matrix_view_array (Iij, 3, 3);
    //gsl_vector *eval = gsl_vector_alloc (3);
    //gsl_matrix *evec = gsl_matrix_alloc (3, 3);
    //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    //gsl_eigen_symmv (&I.matrix, eval, evec, w);
    //gsl_eigen_symmv_free (w);
    //for (l = 0; l < 3; l++){
    //  for (k = 0; k < 3; k++){
    //    mevac[l*3+k] = gsl_matrix_get (evec, l, k);
    //  }
    //}
    //for (l = 0; l < 3; l++){
    //  evals[l] = SQRT(gsl_vector_get (eval, l));
    //}
    //axisx = 3.0*evals[0]/(evals[0]+evals[1]+evals[2]);
    //axisy = 3.0*evals[1]/(evals[0]+evals[1]+evals[2]);
    //axisz = 3.0*evals[2]/(evals[0]+evals[1]+evals[2]);
    //gsl_matrix_view M = gsl_matrix_view_array (mevac, 3, 3);
    //gsl_vector_free (eval);
    //gsl_matrix_free (evec);

    mergeSortDouble( dist0c, id0c, 0, npinh-1);
    M1 = log(1.+cm) - cm/(1.+cm);
    for (i=0; i< npinh; i++){
      j = *( id0c + i);
      rk = (double) i/ npinh;
      x = -exp(-rk*M1-1);
      y = gsl_sf_lambert_W0(x);//lambert function
      gsl_sf_result result;
      status = gsl_sf_lambert_W0_e(x, &result);
      if (status!=GSL_SUCCESS){
        printf("error in Lambert function\n");
        exit(42);
      }
      q = -1./cm * ( 1. + 1./y );
      q = q * Rvir;
      q*=2.0;
      if (q<0.){
        printf("q should not be negative");
        exit(42);
      }
      r = *(dist0c + i);
      if (r>0.0){
	Dx = *(pos + j ) - poshcx;
	Dy = *( pos + j + nothalonum ) - poshcy;
	Dz = *( pos + j + 2*nothalonum )- poshcz;
        Dx = diffpbc(Dx,boxsize);
        Dy = diffpbc(Dy,boxsize);
        Dz = diffpbc(Dz,boxsize);
        phi = ATAN( Dy , Dx );
        theta = ACOS( Dz/r );
        xyz[0] = q*SIN(theta)*COS(phi);
        xyz[1] = q*SIN(theta)*SIN(phi);
        xyz[2] = q*COS(theta);
	//gsl_vector_view XYZ = gsl_vector_view_array (xyz, 3);
	//gsl_vector_view XYZ2 = gsl_vector_view_array (xyz2, 3);
	//gsl_blas_dgemv ( CblasNoTrans, 1., &M.matrix, &XYZ.vector, 0., &XYZ2.vector);
	//xyz2[0] *= axisx;
	//xyz2[1] *= axisy;
	//xyz2[2] *= axisz;
	//gsl_blas_dgemv ( CblasTrans, 1., &M.matrix, &XYZ2.vector, 0., &XYZ.vector);
        *(pos + j) = pbc(poshcx + xyz[0],boxsize);
        *(pos + j + nothalonum) = pbc(poshcy + xyz[1],boxsize);
        *(pos + j + 2*nothalonum) = pbc(poshcz + xyz[2],boxsize);
      }
    }
    free(id0c);
    free(dist0c);
  }

  /* write to binary */
  sprintf(saveToPath_halocatalogue,"%s_hc.dat", path);
  FILE *hc;
  if ( (hc=fopen(saveToPath_halocatalogue,"wb"))==NULL){
    printf("could not open the halocatalogue file\n");
    exit(0);
  }
  printf("written halo catalogue to file %s_hc.dat\n",path);

  /* number of halos */
  fwrite( &nh, sizeof(int) , 1 , hc );

  for (h=0;h<nh;h++){
    fwrite( (*phalos)[h], sizeof(Halo) + sizeof(int) * (*phalos)[h]->npinh , 1, hc );
  }
  fclose(hc);

  sprintf(saveToPath_halonum,"%s_hnum.dat", path);
  FILE *hc2;
  if ( (hc2=fopen(saveToPath_halonum,"wb"))==NULL){
    printf("could not open the halonum file\n");
    exit(0);
  }

  fwrite( halonum , sizeof(int)*nothalonum , 1 , hc2 );
  printf("written halonum to file %s_hnum.dat\n", path );
  fclose(hc2);

  /* clean up */
  free(numinh);
  free(haloID);
  free(halonum);
  for (h=0;h<nh;h++){
    free(halos[h]);
  }
  return nh;
}
