#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/getneighborlist.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define PI 3.14159265358979323846264338327950288419716939937510
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
#define EC 2.718281828459
#define TEMPPOW 2.0
#define TEMPPOW2 2.0
#define FAC 1.5
#define REPFAC 0.01
#define RANDOMFAC 1
#define PROB 0.2
#define SCALE 0.1
//#define R_SPHERE_SQ 90000











double getrfromlist(double *plist,double *rlist,int prnum,double targetp)
 {
 int i,j,id;
 double targetvalue;

 if(targetp<plist[0])
  targetvalue = rlist[prnum-1];	//if p is too large, two point should be very near, thus ingore the attraction
 else if(targetp > plist[prnum-1])
  targetvalue = rlist[prnum-1];
 else
  {
  id = (targetp-plist[0])/(plist[1]-plist[0]);
  if(id >= prnum)
   id = prnum-1;
  targetvalue = rlist[id]+(rlist[id+1]-rlist[id])*(targetp-plist[id])/(plist[id+1]-plist[id]);
  }

 return(targetvalue);
 }








double randomdouble()
 {
 return rand()/(double)RAND_MAX;
 }


// a normal-distributed rand generator. From Fred & Yutong's code
// using Box-Muller transform (polar form)
double random_norm(double mean, double sigma)
 {
    double r = 2.0;
    double U1 = 0;
    double U2 = 0;

    while (r >= 1.0) {

        U1 = 2.0 * randomdouble() - 1.0;  // [-1, +1]
        U2 = 2.0 * randomdouble() - 1.0;

        r = U1*U1 + U2*U2;

    }

    return mean + sigma * U1 * sqrt(-2.0 * log(r)/ r);

 }







void gen3Drandomvector(double *vector)
 {
 int i,j;
 double theta,phi;

 theta = rand()/(double)RAND_MAX*2-1;
 phi = rand()/(double)RAND_MAX*2*PI;

 theta = acos(theta);

 vector[0] = sin(theta)*cos(phi);
 vector[1] = sin(theta)*sin(phi);
 vector[2] = cos(theta);
 }




//double getdistsq(double *coor1,double *coor2,int dim)
// {
// int i;
// double distsq;
// distsq = 0;
//
// for(i=0;i<dim;i++)
//  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);
//
// return(distsq);
// }
//





int getacceptbool(double **coor,int targetid,int dim,double radius)
 {
 int acceptbool;
 int i,j,k;
 double distsq;

 acceptbool = 1;
 for(i=0;(i<targetid)&&(acceptbool==1);i++)
  {
  distsq = getdistsq(coor[i],coor[targetid],dim);
  if(distsq < 4*radius*radius)
   acceptbool = 0;
  }

 return(acceptbool);

 }






double getdiff(double data1,double data2)
 {
 double diff;
 diff = data1-data2;
 if(diff < 0)
  diff = -diff;

 return(diff);
 }




double getnorm(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 norm = sqrt(norm);
 return(norm);
 
 }




double getnormsq(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 return(norm);
 
 }





void normalizevector(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 if(norm < 0)
  {
  printf("Errer when normaling vectors\n");
  exit(0);
  }

 norm = sqrt(norm);
 for(i=0;i<dim;i++)
  vector[i] /= norm;

 }







void genvectormatrix(double **coor,int statenum,double ***vectorlist)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(i<j)
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = coor[j][k]-coor[i][k];
    normalizevector(vectorlist[i][j],DIM);
    }
   else
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = vectorlist[j][i][k];
   }
 for(i=0;i<statenum;i++)
  for(k=0;k<DIM;k++)
   vectorlist[i][i][k] = 0;
 }








void genvectordistmatrix(double **coor,int statenum,double ***vectorlist,double **distmatrix)
 {
 int i,j,k;


 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(i<j)
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = coor[j][k]-coor[i][k];
    distmatrix[i][j] = getnorm(vectorlist[i][j],DIM);
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] /= distmatrix[i][j];
//    normalizevector(vectorlist[i][j],DIM);
    }
   else
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = -vectorlist[j][i][k];
    distmatrix[i][j] = distmatrix[j][i];
    }
   }

 for(i=0;i<statenum;i++)
  {
  distmatrix[i][i] = 0;
  for(k=0;k<DIM;k++)
   vectorlist[i][i][k] = 0;
  }
 }





void genvectordistmatrix_block(double **coor,int statenum,double ***vectormatrix,double **distmatrix,double **eqrmatrix)
 {
 int i,j,k;
 double distmax;
 distmax = 90000;


 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(eqrmatrix[i][j] < distmax)
    {
    for(k=0;k<DIM;k++)
     vectormatrix[i][j][k] = coor[j][k]-coor[i][k];
    distmatrix[i][j] = getnorm(vectormatrix[i][j],DIM);
    for(k=0;k<DIM;k++)
     vectormatrix[i][j][k] /= distmatrix[i][j];
//    normalizevector(vectorlist[i][j],DIM);

    for(k=0;k<DIM;k++)
     vectormatrix[j][i][k] = -vectormatrix[i][j][k];
    distmatrix[j][i] = distmatrix[i][j];
    }
   }

// for(i=0;i<statenum;i++)
//  {
//  distmatrix[i][i] = 0;
//  for(k=0;k<DIM;k++)
//   vectorlist[i][i][k] = 0;
//  }
 }






void genvectordistlist(double **coor,int statenum,double **vectorlist,int *id1list,int *id2list,double *distlist,int nonzeronum,double **vectorlist_neig,double *distlist_neig)
 {
 int i,j,k,id1,id2;
 double dist,tempvector[DIM];

 for(i=0;i<nonzeronum;i++)
  {
  id1 = id1list[i];
  id2 = id2list[i];
  for(k=0;k<DIM;k++)
   vectorlist[i][k] = coor[id2][k]-coor[id1][k];
  distlist[i] = getnorm(vectorlist[i],DIM);
  for(k=0;k<DIM;k++)
   vectorlist[i][k] /= distlist[i];
  }


 for(i=0;i<statenum-1;i++)
  {
  for(k=0;k<DIM;k++)
   vectorlist_neig[i][k] = coor[i+1][k]-coor[i][k];
  dist = getnorm(vectorlist_neig[i],DIM);
  for(k=0;k<DIM;k++)
   vectorlist_neig[i][k] /=dist;
  distlist_neig[i] = dist;
  }

 }












void modifytwopoint(double **coor,int statenum,int id1,int id2,double radius)
 {
 int i,j,k;
 double vector[DIM],dist,center[DIM];

 dist=0;
 for(i=0;i<DIM;i++)
  {
  vector[i] = coor[id2][i]-coor[id1][i];
  dist += vector[i]*vector[i];
  center[i] = (coor[id1][i]+coor[id2][i])/2;
  }

 dist=sqrt(dist);
 for(i=0;i<DIM;i++)
  vector[i] /= dist;

 for(i=0;i<DIM;i++)
  {
  coor[id1][i] = center[i]-vector[i]*radius;
  coor[id2][i] = center[i]+vector[i]*radius;
  }


 }








int getcollisionbool(double **coor,int statenum,double radius)
 {
 int i,j,k,collisionbool;
 double distsq,cutoff;

 cutoff = 4*radius*radius;

 collisionbool = 0;
 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < cutoff)
    {
    collisionbool = 1;
    modifytwopoint(coor,statenum,i,j,radius);
    }
   }
 return(collisionbool);
 }









void getforcematrix(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double **alphamatrix,double bondfactor,double bondlen)
 {
 int i,j,k;
 double forcesq,temp;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    temp = (-alphamatrix[i][k]/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
    forcelist[i][j] += temp;
    forcelist[k][j] -= temp;
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }









void getforce(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double alpha,double bondfactor,double bondlen)
 {
 int i,j,k;
 double forcesq,temp;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
    forcelist[i][j] += temp;
    forcelist[k][j] -= temp;
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }




// getforcematrix_harm(vectorlist,distmatrix,coefmatrix,statenum,forcelist,bondfactor,bondlen,eqrmatrix,drmatrix,slopematrix);	//bonded energy=kb(r-r0)^2, nonbonded energy = alpha/r+kr

void getforcematrix_harm(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **drmatrix,double **slopematrix)
 {
 int i,j,k;
 double forcesq,temp;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    if(distmatrix[i][k] > rmatrix[i][k]+drmatrix[i][k])
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }






//void getforcematrix_harm(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **drmatrix,double **slopematrix)
void getforcematrix_Morse(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix)
 {
 int i,j,k;
 double forcesq,temp;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    if(distmatrix[i][k] > rmatrix[i][k]*1.2)
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
//     temp = 2*De*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*vectorlist[i][k][j]/rmatrix[i][k];	//a*De*(1-e^(2(r-re)/r2))
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*4/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
     temp = 2*Dematrix[i][k]*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*2/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*1/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=1/re
//     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }







void getforcematrix_Morse_noneighbour(double ***vectormatrix,double **distmatrix,int statenum,double **forcematrix,double **eqrmatrix,double **coefmatrix,double **slopematrix,double *powlist,double *elist,int elen)
 {
 int i,j,k,id1,id2;
 double forcesq,temp,temp2,temp3,distmax;

 distmax = 90000;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcematrix[i][j] = 0;


 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   if(eqrmatrix[i][j] < distmax)
    {
    if(distmatrix[i][j] > eqrmatrix[i][j]*1.5)
     {
     temp = slopematrix[i][j];
     for(k=0;k<DIM;k++)
      {
      forcematrix[i][k] += temp*vectormatrix[i][j][k];
      forcematrix[j][k] += temp*vectormatrix[j][i][k];
      }
     }
    else if(distmatrix[i][j] < eqrmatrix[i][j]*0.8)
     {
     temp = slopematrix[i][j]*50*(eqrmatrix[i][j]*eqrmatrix[i][j]/(distmatrix[i][j]*distmatrix[i][j])-1);
     for(k=0;k<DIM;k++)
      {
      forcematrix[i][k] -= temp*vectormatrix[i][j][k];
      forcematrix[j][k] -= temp*vectormatrix[j][i][k];
      }
     }
    else
     {
     temp2=-8*(distmatrix[i][j]/eqrmatrix[i][j]-1);
     temp3 = getrfromlist(powlist,elist,elen,temp2);
     temp = 2*coefmatrix[i][j]*(1-temp3)*temp3*8/eqrmatrix[i][j];
  
  
     for(k=0;k<DIM;k++)   
      {
  //    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
      forcematrix[i][k] += temp*vectormatrix[i][j][k];
      forcematrix[j][k] += temp*vectormatrix[j][i][k];
      }
     }
    }
   }





 }



//void getforcematrix_Morse_noneighbour_all(int blocknum,double **blockeqrmatrix,double **blockcoefmatrix,double **blockslopematrix,double *powlist,double *elist,int elen,double **forcelist,int statenum,int *assign,double **coor,double distmax)
double getforcematrix_Morse_noneighbour_all(int blocknum,double **blockeqrmatrix,double **blockcoefmatrix,double **blockslopematrix,double *powlist,double *elist,int elen,double **forcelist,int statenum,int *assign,double **coor,double distmax,double factor)
 {
 int i,j,k,id1,id2;
 double dist,*vector,temp,temp2,temp3,Fmax;

 vector = doublearray(DIM);
 distmax *= 0.9;

 Fmax = 0;
 for(i=0;i<statenum;i++)
  {
  id1 = assign[i];
  for(j=i;(j<statenum)&&(id1!=-1);j++)	//only consider half of the block matrix
   {
   id2 = assign[j];
   if((id2!=-1)&&(blockeqrmatrix[id1][id2] < distmax))
    {
    for(k=0;k<DIM;k++)
     vector[k] = coor[j][k]-coor[i][k];
    dist = getnorm(vector,DIM);
    for(k=0;k<DIM;k++)
     vector[k] /= dist;

    if(dist > blockeqrmatrix[id1][id2]*1.0)
     {
//     temp = slopelist[i]*distlist[i]/eqrlist[i];
//     temp = 0.1*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
//     temp = FAC*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
     temp = factor*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += temp*vector[k];
      forcelist[j][k] -= temp*vector[k];
      }
     }
    else if((dist < blockeqrmatrix[id1][id2]*1.0)&&(dist > blockeqrmatrix[id1][id2]*0.0))
     {
//     temp = blockslopematrix[id1][id2]*40*dist/blockeqrmatrix[id1][id2];
//     temp = blockslopematrix[id1][id2]*50;
     temp = blockslopematrix[id1][id2]*50*(blockeqrmatrix[id1][id2]/(dist)-1);
//     temp = blockslopematrix[id1][id2]*50*(blockeqrmatrix[id1][id2]*blockeqrmatrix[id1][id2]/(dist*dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] -= temp*vector[k];
      forcelist[j][k] += temp*vector[k];
      }
     }
//    else if((dist < blockeqrmatrix[id1][id2]*1.5)&&(dist > blockeqrmatrix[id1][id2]*0.8))
//     {
//     temp2=-8*(dist/blockeqrmatrix[id1][id2]-1);
////     temp3 = getrfromlist(powlist,elist,elen,temp2);
//     temp3 = pow(EC,temp2);
//     temp = 2*blockcoefmatrix[id1][id2]*(1-temp3)*temp3*8/blockeqrmatrix[id1][id2];
//  
//  
//     for(k=0;k<DIM;k++)   
//      {
//  //    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
//      forcelist[i][k] += temp*vector[k];
//      forcelist[j][k] -= temp*vector[k];
//      }
//     }
  if(temp < 0)
   temp *= -1;
  if(Fmax < temp)
   Fmax = temp;
//   printf("block %d %d %lf\n",i,j,temp);
 
    }
   }
  }


 free(vector);
 
 }





// getforcematrix_Morsefromlist(vectorlist,distmatrix,statenum,forcelist,bondfactor,bondlen,eqrmatrix,coefmatrix,slopematrix,powlist,elist,elen);	//De=kT,a=2/re,dr=re
void getforcematrix_Morsefromlist(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix,double *powlist,double *elist,int elen)
 {
 int i,j,k;
 double forcesq,temp,temp2,temp3;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    if(distmatrix[i][k] > rmatrix[i][k]*3)
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
     temp2=-2*(distmatrix[i][k]/rmatrix[i][k]-1);
     temp3 = getrfromlist(powlist,elist,elen,temp2);
     
//     temp = 2*De*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*vectorlist[i][k][j]/rmatrix[i][k];	//a*De*(1-e^(2(r-re)/r2))
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*4/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
     temp = 2*Dematrix[i][k]*(1-temp3)*temp3*vectorlist[i][k][j]*2/rmatrix[i][k];
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*2/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=2/re
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*1/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=1/re
//     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }





//void getforcematrix_Morsefromlist(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix,double *powlist,double *elist,int elen)

//void getforcelist_Morsefromlist(double **vectorlist,double *distlist,int statenum,int *id1list,int *id2list,int nonzeronum,double **forcelist,double bondfactor,double bondlen,double *eqrlist,double *coeflist,double *slopelist,double *powlist,double *elist,int elen,double **vectorlist_neig,double *distlist_neig)
double getforcelist_Morsefromlist(double **vectorlist,double *distlist,int statenum,int *id1list,int *id2list,int nonzeronum,double **forcelist,double bondfactor,double bondlen,double *eqrlist,double *coeflist,double *slopelist,double *powlist,double *elist,int elen,double **vectorlist_neig,double *distlist_neig,double factor)
 {
 int i,j,k,id1,id2;
 double forcesq,temp,temp2,temp3,Fmax;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

Fmax = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(distlist_neig[i]-bondlen)*vectorlist_neig[i][j];
   forcelist[i+1][j]  -= bondfactor*(distlist_neig[i]-bondlen)*vectorlist_neig[i][j];
   }
  }

//calculate nonbonded force
 for(i=0;i<nonzeronum;i++)
  {
  id1 = id1list[i];
  id2 = id2list[i];
  if(distlist[i] >= eqrlist[i]*1.0)
   {
//    temp = 0.1*slopelist[i]*distlist[i]/eqrlist[i];
//    temp = FAC*slopelist[i]*distlist[i]/eqrlist[i];
    temp = factor*slopelist[i]*(distlist[i]/eqrlist[i]-1);
//    temp = slopelist[i];
//    temp = slopelist[i]*(1+0.3*distlist[i]/eqrlist[i]);
   for(k=0;k<DIM;k++)
    {
//    temp = slopelist[i]*vectorlist[i][k];
    forcelist[id1][k] += temp*vectorlist[i][k];
    forcelist[id2][k] -= temp*vectorlist[i][k];
    }
   }
  else if((distlist[i] < eqrlist[i]*1.0)&&(distlist[i] > eqrlist[i]*0.0))
   {
//   temp = slopelist[i]*40*distlist[i]/eqrlist[i];
//   temp = slopelist[i]*50;
   temp = slopelist[i]*50*(eqrlist[i]/(distlist[i])-1);
//   temp = slopelist[i]*50*(eqrlist[i]*eqrlist[i]/(distlist[i]*distlist[i])-1);
   for(k=0;k<DIM;k++)
    {
//    temp = slopelist[i]*vectorlist[i][k]*eqr;
    forcelist[id1][k] -= temp*vectorlist[i][k];
    forcelist[id2][k] += temp*vectorlist[i][k];
    }
   }
//  else if((distlist[i] < eqrlist[i]*1.5)&&(distlist[i] > eqrlist[i]*0.8))
//   {
//   temp2=-8*(distlist[i]/eqrlist[i]-1);
//   temp3 = getrfromlist(powlist,elist,elen,temp2);
//   temp = 2*coeflist[i]*(1-temp3)*temp3*8/eqrlist[i];
//
//
//   for(k=0;k<DIM;k++)   
//    {
////    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
//    forcelist[id1][k] += temp*vectorlist[i][k];
//    forcelist[id2][k] -= temp*vectorlist[i][k];
//    }
//   }
  if(temp < 0)
   temp *= -1;
  if(Fmax < temp)
   Fmax = temp;
//  printf("nonblock %d %d %lf\n",id1,id2,temp);
  }

 return(Fmax);

 }


















void genvel(double **velvector,int statenum,int dim,double sigma)
 {
 int i,j;
 for(i=0;i<statenum;i++)
  for(j=0;j<dim;j++)
   velvector[i][j] = random_norm(0,sqrt(sigma));
 }






void genvel_mass(double **velvector,int statenum,int dim,double sigma,double *mass)
 {
 int i,j;
 for(i=0;i<statenum;i++)
  for(j=0;j<dim;j++)
   velvector[i][j] = random_norm(0,sqrt(sigma/mass[i]));
 }











void getfinalvel_mass(double **velvector,int statenum,double **forcelist,double dt,double **velvector2,double *mass)
 {
 int i,j,k;

 for(i=0;i<statenum;i++) 
  for(j=0;j<DIM;j++)
   velvector2[i][j] = velvector[i][j]+forcelist[i][j]*dt/mass[i];
 }







void getfinalvel(double **velvector,int statenum,double **forcelist,double dt,double **velvector2)
 {
 int i,j,k;

 for(i=0;i<statenum;i++) 
  for(j=0;j<DIM;j++)
   velvector2[i][j] = velvector[i][j]+forcelist[i][j]*dt;
 }








void updatecoor(double **coor,double **velvector,double **velvector2,int statenum,double dt)
 {
 int i,j;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] += (velvector[i][j]+velvector2[i][j])*dt/2;

 }



void copyvel(double **velvector2,double **velvector,int statenum,int dim)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  for(j=0;j<dim;j++)
  velvector[i][j] = velvector2[i][j];

 }







void genblockcenter(double **coor,int statenum,int *assign,double **blockcentercoor,int blocknum,int *blocknumlist)
 {
 int i,j,k;

 for(i=0;i<blocknum;i++)
  for(j=0;j<DIM;j++)
   blockcentercoor[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  if(assign[i] != -1)
   {
   for(j=0;j<DIM;j++)
    blockcentercoor[assign[i]][j] += coor[i][j];
   }
  }

 for(i=0;i<blocknum;i++)
  {
  if(blocknumlist[i] != 0)
   {
   for(j=0;j<DIM;j++)
    blockcentercoor[i][j] /= blocknumlist[i];
   }
  }
 }





void assignblockforce(double **blockforcematrix,int blocknum,double **forcelist,int statenum,int *assign)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  {
  if(assign[i] != -1)
   {
   for(j=0;j<DIM;j++)
    forcelist[i][j] += blockforcematrix[assign[i]][j];
   }
  }


 }





void getinitialpairing(double **coor,int chrlen,int *chrtype,double **pointcoor,int pointnum,int *pointtype,int *chrcontactid,int *pointcontactid,double nearcutoff,double contactprob)
 {
 int i,j,k,endbool,nearbool;
 double distsq,nearcutoffsq;

 nearcutoffsq = nearcutoff*nearcutoff;
 for(i=0;i<chrlen;i++)
  {
  endbool = 0;
  for(j=0;(j<pointnum)&&(endbool == 0);j++)
   {
   if((pointtype[j]==chrtype[i])&&(pointcontactid[j] == -1))
    {
    distsq = 0;
    nearbool = 1;
    for(k=0;(k<DIM)&&(nearbool == 1);k++)
     {
     distsq += (coor[i][k]-pointcoor[j][k])*(coor[i][k]-pointcoor[j][k]);
     if(distsq > nearcutoffsq)
      nearbool = 0;
     }

    if(nearbool == 1)
     {
     if(randomdouble()<contactprob)
      {
      chrcontactid[i] = j;
      pointcontactid[j] = i;
      endbool = 1;
      }
     }
    }
   }
  }


 }





void updatepairing(double **coor,int chrlen,int *chrtype,double **pointcoor,int pointnum,int *pointtype,int *chrcontactid,int *pointcontactid,double nearcutoff,double contactprob,double sepprob)
 {
 int i,j,k,endbool,nearbool;
 double distsq,nearcutoffsq;

 nearcutoffsq = nearcutoff*nearcutoff;


 for(i=0;i<chrlen;i++)
  {
  if(chrcontactid[i] != -1)
   {
   if(randomdouble()<sepprob)
    {
    pointcontactid[chrcontactid[i]] = -1;
    chrcontactid[i] = -1;
    }
   }

  if(chrcontactid[i] == -1)
   {
   endbool = 0;
   for(j=0;(j<pointnum)&&(endbool == 0);j++)
    {
    if((pointtype[j]==chrtype[i])&&(pointcontactid[j] == -1))
     {
     distsq = 0;
     nearbool = 1;
     for(k=0;(k<DIM)&&(nearbool == 1);k++)
      {
      distsq += (coor[i][k]-pointcoor[j][k])*(coor[i][k]-pointcoor[j][k]);
      if(distsq > nearcutoffsq)
       nearbool = 0;
      }
 
     if(nearbool == 1)
      {
      if(randomdouble()<contactprob)
       {
       chrcontactid[i] = j;
       pointcontactid[j] = i;
       endbool = 1;
       }
      }
     }
    }
   }
  }


 }






double getforce_v2(double **coor,int chrlen,double **pointcoor,int pointnum,int *chrcontactid,int *pointcontactid,double nearcutoff,double kb,double **chrforce,double **pointforce,double beadradius,double bondlen,double sphereradius,double pointmass)
 {
 int i,j,k,l,tempid;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax;

 coef = kb*bondlen*bondlen/(nearcutoff*nearcutoff);	//automatrix scale with factor of 0.25
 tempvector = doublearray(DIM);
 tempcutoff = beadradius*beadradius;

 for(i=0;i<chrlen;i++)
  for(j=0;j<DIM;j++)
   chrforce[i][j] = 0;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   pointforce[i][j] = 0;

 //getsphere restriction

 for(i=0;i<chrlen;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
   for(j=0;j<DIM;j++)
    chrforce[i][j] = -coor[i][j]/sphereradius;
   }
  }

 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(pointcoor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
   for(j=0;j<DIM;j++)
    pointforce[i][j] = -pointcoor[i][j]/sphereradius;
   }
  }

 //get bonded interaction
 for(i=0;i<chrlen-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);
  for(j=0;j<DIM;j++)
   {
   chrforce[i][j] += force*tempvector[j]/dist;
   chrforce[i+1][j] -= force*tempvector[j]/dist;
   }
  }

 for(i=0;i<chrlen-2;i++)
  {
  for(j=i+2;j<chrlen;j++)
   {
   dist = getdistsq(coor[i],coor[j],DIM);
   if(dist < beadradius)
    {
    force = 1*kb*bondlen*bondlen/dist;
    for(k=0;k<DIM;k++)
     {
     chrforce[i][k] -= (coor[j][k]-coor[i][k])*force/dist;
     chrforce[j][k] += (coor[j][k]-coor[i][k])*force/dist;
     }
    }
   }
  }


//get nonbonded interaction
 for(i=0;i<chrlen;i++)
  {
  if(chrcontactid[i] != -1)
   {
   tempid = chrcontactid[i];
   dist = getdistsq(coor[i],pointcoor[tempid],DIM);
   dist = sqrt(dist);
   for(j=0;j<DIM;j++)
    tempvector[j] = pointcoor[tempid][j]-coor[i][j];
 
   force = 2*coef*(dist-nearcutoff/2);
   for(j=0;j<DIM;j++)
    {
    chrforce[i][j] += force*tempvector[j]/dist;
    pointforce[tempid][j] -= force*tempvector[j]/dist;
    }
   }
  }

 Fmax = 0;
 for(i=0;i<chrlen;i++)
  {
  force = getnormsq(chrforce[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(pointforce[i],DIM);
  if(force/pointmass > Fmax)
   Fmax = force/pointmass;
  }

 free(tempvector);
 return(Fmax);
 }






void updatecoor_v2(double **coor,int chrlen,double **pointcoor,int pointnum,double **chrvel,double **pointvel,double **chrforce,double **pointforce,double dt,double pointmass)
 {
 int i,j,k,l;
 double tempcoor;

 for(i=0;i<chrlen;i++)
  {
  for(j=0;j<DIM;j++)
   {
   coor[i][j] += chrvel[i][j]*dt+0.5*chrforce[i][j]*dt*dt;
   }
  }

 for(i=0;i<pointnum;i++)
  {
  for(j=0;j<DIM;j++)
   pointcoor[i][j] += pointvel[i][j]*dt+0.5*pointforce[i][j]*dt*dt/pointmass;
  }


 }




double getgyr(double **coor,int chrlen)
 {
 int i,j;
 double *center,gyrsq;

 center = doublearray(DIM);

 for(i=0;i<DIM;i++)
  center[i] = 0;

 for(i=0;i<chrlen;i++)
  {
  for(j=0;j<DIM;j++)
   center[j] += coor[i][j];
  }

 for(j=0;j<DIM;j++)
  center[j] /= chrlen;

 gyrsq = 0;
 for(i=0;i<chrlen;i++)
  {
  for(j=0;j<DIM;j++)
   gyrsq += (coor[i][j]-center[j])*(coor[i][j]-center[j]);
  }

 gyrsq = sqrt(gyrsq/chrlen);


 }






void getfeatureassign(int *featureassign,int featurenum,int assignvalue)
 {
 int i,j,k,remainvalue;

 for(i=0;i<featurenum;i++)
  {
  featureassign[i] = assignvalue%2;
  assignvalue /= 2;
  }
 }





void getfeatureassign_randomtype(int *featureassign,int featurenum,int assignvalue)
 {
 int i,j,k,remainvalue;

 for(i=0;i<featurenum;i++)
  {
  featureassign[i] = assignvalue%2;
  assignvalue /= 2;
  }

 featureassign[featurenum-1] = 1;
 for(i=0;i<featurenum-1;i++)
  if(featureassign[i] == 1)
   featureassign[featurenum-1] = 0;
 }







void gettargetlist(int **featureassign,int pointnum,int featurenum,int **targetlist,int *targetlen)
 {
 int i,j,k;
 for(i=0;i<featurenum;i++)
  targetlen[i] = 0;

 for(j=0;j<featurenum;j++)
  {
  for(i=0;i<pointnum;i++)
   {
   if(featureassign[i][j] == 1)
    {
    targetlist[j][targetlen[j]] = i;
    targetlen[j] ++;
    }
   }
  }

 }



void getpairlist(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,double distcutoff)
 {
 int i,j,k,l,id1,id2;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  pairlen[i] = 0;

 for(i=0;i<featurenum;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+3))
     {
     pairlist[id1][pairlen[id1]] = id2;
     pairlen[id1] ++;
     }
    }
   }
  }


 }
 




void getpairlist_random(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,double distcutoff)
 {
 int i,j,k,l,id1,id2;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  randompairlen[i] = 0;
  }


 for(i=0;i<featurenum-1;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+1))
     {
     pairlist[id1][pairlen[id1]] = id2;
     pairlen[id1] ++;
     }
    }
   }
  }


 for(i=featurenum-1;i<=featurenum-1;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+1))
     {
     randompairlist[id1][pairlen[id1]] = id2;
     randompairlen[id1] ++;
     }
    }
   }
  }


 }
 




void getpairlist_random_allagg(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,double distcutoff,int *assign,int **featureassign)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  randompairlen[i] = 0;
  }



 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;(j<pointnum)&&(assign[i] != 0);j++)
   {
   if(assign[j] != 0)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
    if(distsq < distcutoffsq)
     {
     pairlist[i][pairlen[i]] = j;
     pairlen[i] ++;
     }
    }
   }
  }




// for(i=0;i<featurenum-1;i++)
//  {
//  for(j=0;j<targetlen[i];j++)
//   {
//   id1 = targetlist[i][j];
//   for(k=j+1;k<targetlen[i];k++)
//    {
//    id2 = targetlist[i][k];
//    distsq = getdistsq(coor[id1],coor[id2],DIM);
//    if((distsq < distcutoffsq)&&(id2 > id1+1))
//     {
//     pairlist[id1][pairlen[id1]] = id2;
//     pairlen[id1] ++;
//     }
//    }
//   }
//  }


 for(i=featurenum-1;i<=featurenum-1;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+1))
     {
     randompairlist[id1][pairlen[id1]] = id2;
     randompairlen[id1] ++;
     }
    }
   }
  }


 }
 




void getpairlist_random_allagg_prob(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,double distcutoff,int *assign,int **featureassign)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  randompairlen[i] = 0;
  }



 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;(j<pointnum)&&(assign[i] != 0);j++)
   {
   if(assign[j] != 0)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
    if(distsq < distcutoffsq)
     {
     endbool = 0;
     for(k=0;(k<featurenum-1)&&(endbool==0);k++)
      {
      for(l=0;(l<featurenum-1)&&(endbool==0)&&(featureassign[i][k] == 1);l++)
       {
       if((randomdouble() < PROB)&&(featureassign[j][l]==1))
        {
        pairlist[i][pairlen[i]] = j;
        pairlen[i] ++;
        endbool = 1;
        }
       }
      }
//     pairlist[i][pairlen[i]] = j;
//     pairlen[i] ++;
     }
    }
   }
  }




// for(i=0;i<featurenum-1;i++)
//  {
//  for(j=0;j<targetlen[i];j++)
//   {
//   id1 = targetlist[i][j];
//   for(k=j+1;k<targetlen[i];k++)
//    {
//    id2 = targetlist[i][k];
//    distsq = getdistsq(coor[id1],coor[id2],DIM);
//    if((distsq < distcutoffsq)&&(id2 > id1+1))
//     {
//     pairlist[id1][pairlen[id1]] = id2;
//     pairlen[id1] ++;
//     }
//    }
//   }
//  }


 for(i=featurenum-1;i<=featurenum-1;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+1))
     {
     randompairlist[id1][pairlen[id1]] = id2;
     randompairlen[id1] ++;
     }
    }
   }
  }


 }
 




void getpairlist_random_allagg_prob_time(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,double distcutoff,int *assign,int **featureassign,int **contacttime)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  randompairlen[i] = 0;
  for(j=0;j<pointnum;j++)
   contacttime[i][j] = 0;
  }



 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;(j<pointnum)&&(assign[i] != 0);j++)
   {
   if(assign[j] != 0)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
    if(distsq < distcutoffsq)
     {
     endbool = 0;
     for(k=0;(k<featurenum-1);k++)
      {
      for(l=0;(l<featurenum-1)&&(featureassign[i][k] == 1);l++)
       {
       if((randomdouble() < PROB)&&(featureassign[j][l]==1))
        {
        if(endbool == 0)
         {
         pairlen[i]++;
         pairlist[i][pairlen[i]-1] = j;
         endbool = 1;
         contacttime[i][pairlen[i]-1] ++;
         }
        else
         {
         pairlist[i][pairlen[i]-1] = j;
         contacttime[i][pairlen[i]-1] ++;
         }
//        pairlist[i][pairlen[i]] = j;
//        pairlen[i] ++;
//        endbool = 1;
        }
       }
      }
//     pairlist[i][pairlen[i]] = j;
//     pairlen[i] ++;
     }
    }
   }
  }




// for(i=0;i<featurenum-1;i++)
//  {
//  for(j=0;j<targetlen[i];j++)
//   {
//   id1 = targetlist[i][j];
//   for(k=j+1;k<targetlen[i];k++)
//    {
//    id2 = targetlist[i][k];
//    distsq = getdistsq(coor[id1],coor[id2],DIM);
//    if((distsq < distcutoffsq)&&(id2 > id1+1))
//     {
//     pairlist[id1][pairlen[id1]] = id2;
//     pairlen[id1] ++;
//     }
//    }
//   }
//  }


 for(i=featurenum-1;i<=featurenum-1;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+3)&&(randomdouble()<PROB))
     {
     randompairlist[id1][pairlen[id1]] = id2;
     randompairlen[id1] ++;
     }
    }
   }
  }


 }
 








void getpairlist_random_allagg_prob_time_alllist_v2(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,double distcutoff, int *assign,int **featureassign,int **allpairlist,int *allpairlen,int *tempassign)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;


 for(i=0;i<pointnum;i++)
  {
// printf("test 1\n");
  free(pairlist[i]);
// printf("test 2\n");
  free(allpairlist[i]);
  }


 getneighborlist_list_all_split_3_allpair(coor,pointnum,tempassign,allpairlist,allpairlen,distcutoff);


 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  pairlist[i] = intarray(allpairlen[i]);
//  for(j=0;j<allpairlen[i];j++)
//   {
//   id1 = allpairlist[i][j];
//   endbool = 0;
//   for(k=0;(k<featurenum)&&(endbool == 0);k++)
//    if((featureassign[i][k] == 1)&&(featureassign[id1][k] == 1))
//     {
//     pairlen[i]++;
//     pairlist[i][pairlen[i]-1] = id1;
//     endbool = 1;
//     }
//   }
  }

// for(i=0;i<pointnum;i++)
//  {
//// printf("test 1\n");
////  free(pairlist[i]);
//// printf("test 2\n");
//  free(allpairlist[i]);
//  }
//
// getneighborlist_list_all_split_3_allpair(coor,pointnum,tempassign,allpairlist,allpairlen,distcutoff);



// distcutoffsq = distcutoff*distcutoff;
//
// for(i=0;i<pointnum;i++)
//  {
//  pairlen[i] = 0;
////  randompairlen[i] = 0;
//  allpairlen[i] = 0;
////  for(j=0;j<pointnum;j++)
////   contacttime[i][j] = 0;
//  }
//
//
//
//
// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;(j<pointnum);j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < distcutoffsq)
//    {
//    allpairlist[i][allpairlen[i]] = j;
//    allpairlen[i] ++;
//    }
//
//   if((assign[j] != 0)&&(assign[i] != 0))
//    {
////    distsq = getdistsq(coor[i],coor[j],DIM);
//    if(distsq < distcutoffsq)
//     {
//     endbool = 0;
//     for(k=0;(k<featurenum)&&(endbool == 0);k++)
//      {
//      if((featureassign[j][k] == 1)&&(featureassign[i][k] == 1))
//       {
//       if(endbool == 0)
//        {
//        pairlen[i]++;
//        pairlist[i][pairlen[i]-1] = j;
//        endbool = 1;
////        contacttime[i][pairlen[i]-1] ++;
//        }
////       else
////        {
////        pairlist[i][pairlen[i]-1] = j;
////        contacttime[i][pairlen[i]-1] ++;
////        }
//       }
////      for(l=0;(l<featurenum)&&(featureassign[i][k] == 1);l++)
////       {
////       if((randomdouble() < PROB)&&(featureassign[j][l]==1))
////        {
////        if(endbool == 0)
////         {
////         pairlen[i]++;
////         pairlist[i][pairlen[i]-1] = j;
////         endbool = 1;
////         contacttime[i][pairlen[i]-1] ++;
////         }
////        else
////         {
////         pairlist[i][pairlen[i]-1] = j;
////         contacttime[i][pairlen[i]-1] ++;
////         }
//////        pairlist[i][pairlen[i]] = j;
//////        pairlen[i] ++;
//////        endbool = 1;
////        }
////       }
//      }
////     pairlist[i][pairlen[i]] = j;
////     pairlen[i] ++;
//     }
//    }
//   }
//  }
//
//
//
//
//
//
//// for(i=featurenum-1;i<=featurenum-1;i++)
////  {
////  for(j=0;j<targetlen[i];j++)
////   {
////   id1 = targetlist[i][j];
////   for(k=j+1;k<targetlen[i];k++)
////    {
////    id2 = targetlist[i][k];
////    distsq = getdistsq(coor[id1],coor[id2],DIM);
////    if((distsq < distcutoffsq)&&(id2 > id1+3)&&(randomdouble()<PROB))
////     {
////     randompairlist[id1][pairlen[id1]] = id2;
////     randompairlen[id1] ++;
////     }
////    }
////   }
////  }
////

 }




 



void getpairlist_random_allagg_prob_time_alllist_v3(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,double distcutoff,double sphereradius,int *assign,int **featureassign,int **allpairlist,int *allpairlen,int *tempassign)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;


 for(i=0;i<pointnum;i++)
  {
// printf("test 1\n");
  free(pairlist[i]);
// printf("test 2\n");
  free(allpairlist[i]);
  }


 getneighborlist_list_all_split_3_allpair(coor,pointnum,tempassign,allpairlist,allpairlen,sphereradius*0.2);


 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  pairlist[i] = intarray(allpairlen[i]);
  for(j=0;j<allpairlen[i];j++)
   {
   id1 = allpairlist[i][j];
   endbool = 0;
   for(k=0;(k<featurenum)&&(endbool == 0);k++)
    if((featureassign[i][k] == 1)&&(featureassign[id1][k] == 1))
     {
     pairlen[i]++;
     pairlist[i][pairlen[i]-1] = id1;
     endbool = 1;
     }
   }
  }

 for(i=0;i<pointnum;i++)
  {
  free(allpairlist[i]);
  }

 getneighborlist_list_all_split_3_allpair(coor,pointnum,tempassign,allpairlist,allpairlen,distcutoff);


 }




 




void getpairlist_random_allagg_prob_time_alllist(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,double distcutoff,int *assign,int **featureassign,int **contacttime,int **allpairlist,int *allpairlen)
 {
 int i,j,k,l,id1,id2,endbool,targetbool;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  {
  pairlen[i] = 0;
  randompairlen[i] = 0;
  allpairlen[i] = 0;
  for(j=0;j<pointnum;j++)
   contacttime[i][j] = 0;
  }




 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;(j<pointnum);j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < distcutoffsq)
    {
    allpairlist[i][allpairlen[i]] = j;
    allpairlen[i] ++;
    }

   if((assign[j] != 0)&&(assign[i] != 0))
    {
//    distsq = getdistsq(coor[i],coor[j],DIM);
    if(distsq < distcutoffsq)
     {
     endbool = 0;
     for(k=0;(k<featurenum)&&(endbool == 0);k++)
      {
      if((featureassign[j][k] == 1)&&(featureassign[i][k] == 1))
       {
       if(endbool == 0)
        {
        pairlen[i]++;
        pairlist[i][pairlen[i]-1] = j;
        endbool = 1;
        contacttime[i][pairlen[i]-1] ++;
        }
//       else
//        {
//        pairlist[i][pairlen[i]-1] = j;
//        contacttime[i][pairlen[i]-1] ++;
//        }
       }
//      for(l=0;(l<featurenum)&&(featureassign[i][k] == 1);l++)
//       {
//       if((randomdouble() < PROB)&&(featureassign[j][l]==1))
//        {
//        if(endbool == 0)
//         {
//         pairlen[i]++;
//         pairlist[i][pairlen[i]-1] = j;
//         endbool = 1;
//         contacttime[i][pairlen[i]-1] ++;
//         }
//        else
//         {
//         pairlist[i][pairlen[i]-1] = j;
//         contacttime[i][pairlen[i]-1] ++;
//         }
////        pairlist[i][pairlen[i]] = j;
////        pairlen[i] ++;
////        endbool = 1;
//        }
//       }
      }
//     pairlist[i][pairlen[i]] = j;
//     pairlen[i] ++;
     }
    }
   }
  }






// for(i=featurenum-1;i<=featurenum-1;i++)
//  {
//  for(j=0;j<targetlen[i];j++)
//   {
//   id1 = targetlist[i][j];
//   for(k=j+1;k<targetlen[i];k++)
//    {
//    id2 = targetlist[i][k];
//    distsq = getdistsq(coor[id1],coor[id2],DIM);
//    if((distsq < distcutoffsq)&&(id2 > id1+3)&&(randomdouble()<PROB))
//     {
//     randompairlist[id1][pairlen[id1]] = id2;
//     randompairlen[id1] ++;
//     }
//    }
//   }
//  }
//

 }




 




void getpairlist_mon(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,double distcutoff)
 {
 int i,j,k,l,id1,id2,randomid;
 double dist,distsq,distcutoffsq,prob;

 prob = 1.0;

 distcutoffsq = distcutoff*distcutoff;

 for(i=0;i<pointnum;i++)
  pairlen[i] = 0;

 for(i=0;i<featurenum;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+3))
     {
     pairlist[id1][pairlen[id1]] = id2;
     pairlen[id1] ++;
     }
    }
   }
  }

 for(i=0;i<pointnum;i++)
  {
  if(pairlen[i] > 0)
   {
   randomid = rand()%pairlen[i];
   pairlist[i][0] =  pairlist[i][randomid];
   pairlen[i] = 1;
   if(randomdouble() > prob)
    pairlen[i] = 0;
   }
  }


 }
 






//void getneighborlist_list_all_split_3(double **coor,int statenum,int *assign,int **neighborlist,int *neighbornum,double fardist)
// {
// int i,j,k,l,m,id1,id2,id3,*templist,nextid,*binnum,**binassign,****binpointid,***binpointnum,dim;
// double dist,distsq,fardistsq,*mincoor,*maxcoor,binsize;
//
// dim =3; 
//
// templist = intarray(statenum);
// 
// mincoor = doublearray(dim);
// maxcoor = doublearray(dim);
// binnum = intarray(dim);
// binassign = intmatrix(statenum,dim);
// binsize = fardist;
//
// for(i=0;i<dim;i++)
//  {
//  mincoor[i] = coor[0][i];
//  maxcoor[i] = coor[0][i];
//  }
//
// for(i=0;i<statenum;i++)
//  {
//  for(j=0;j<dim;j++)
//   {
//   if(mincoor[j] > coor[i][j])
//    mincoor[j] = coor[i][j];
//
//   if(maxcoor[j] < coor[i][j])
//    maxcoor[j] = coor[i][j];
//   }
//  }
//
//
// for(i=0;i<dim;i++)
//  {
//  mincoor[i] -= binsize;
//  maxcoor[i] += binsize;
//  binnum[i] = (maxcoor[i]-mincoor[i])/binsize+1;
//  }
//
// printf("%d %d %d\n",binnum[0],binnum[1],binnum[2]);
// binpointid = intpointmatrixarray(binnum[0],binnum[1],binnum[2]);
// binpointnum = intmatrixarray(binnum[0],binnum[1],binnum[2]);
//
// for(i=0;i<binnum[0];i++)
//  for(j=0;j<binnum[1];j++)
//   for(k=0;k<binnum[2];k++)
//    binpointnum[i][j][k] = 0;
//
// for(i=0;i<statenum;i++)
//  {
//  for(j=0;j<dim;j++)
//   binassign[i][j] = (coor[i][j]-mincoor[j])/binsize;
//  binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]] ++;
//  }
//
// for(i=0;i<binnum[0];i++)
//  for(j=0;j<binnum[1];j++)
//   for(k=0;k<binnum[2];k++)
//    binpointid[i][j][k] = intarray(binpointnum[i][j][k]+1);
//
// for(i=0;i<binnum[0];i++)
//  for(j=0;j<binnum[1];j++)
//   for(k=0;k<binnum[2];k++)
//    binpointnum[i][j][k] = 0;
//
// for(i=0;i<statenum;i++)
//  {
//  for(j=0;j<dim;j++)
//   binassign[i][j] = (coor[i][j]-mincoor[j])/binsize;
//  binpointid[binassign[i][0]][binassign[i][1]][binassign[i][2]][binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]]] = i;
//  binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]] ++;
//  }
// 
// fardistsq = fardist*fardist;
//
//
// //for all pairs
// for(i=0;i<statenum;i++)
//  {
//  neighbornum[i] = 0;
//
//  if(assign[i] != -1)
//   {
//   for(j=-1;j<=1;j++)
//    for(k=-1;k<=1;k++)
//     for(l=-1;l<=1;l++)
//      {
//      id1 = binassign[i][0]+j;
//      id2 = binassign[i][1]+k;
//      id3 = binassign[i][2]+l;
//      if((id1>=0)&&(id1<binnum[0])&&(id2>=0)&&(id2<binnum[1])&&(id3>=0)&&(id3<binnum[2]))
//       {
//       for(m=0;m<binpointnum[id1][id2][id3];m++)
//        {
//        if(assign[binpointid[id1][id2][id3][m]] != -1)
//         {
//         distsq = getdistsq(coor[i],coor[binpointid[id1][id2][id3][m]],dim);
//         if(distsq < fardistsq)
//          {
//          templist[neighbornum[i]] = binpointid[id1][id2][id3][m];
//          neighbornum[i] ++;
//          }
//         }
//        }
//       }
//      }
//   }
//
////  free(neighborlist[i]);
//
////  if(neighbornum[i] != 0)
////   neighborlist[i] = intarray(neighbornum[i]);
//  for(j=0;j<neighbornum[i];j++)
//   neighborlist[i][j] = templist[j];
//  }
//
//// printf("check %d\n",neighborlist[0][0]);
//
// for(i=0;i<binnum[0];i++)
//  {
//  for(j=0;j<binnum[1];j++)
//   {
//   for(k=0;k<binnum[2];k++)
//    free(binpointid[i][j][k]);
//   free(binpointid[i][j]);
//   free(binpointnum[i][j]);
//   }
//  free(binpointid[i]);
//  free(binpointnum[i]);
//  }
// free(binpointid);
// free(binpointnum);
//
// for(i=0;i<statenum;i++)
//  free(binassign[i]);
// free(binassign);
// free(mincoor);
// free(maxcoor);
// free(binnum);
//
// free(templist);
// }












void getpairlist_all(double **coor,int **targetlist,int *targetlen,int featurenum,int pointnum,int **pairlist,int *pairlen,double distcutoff,int **allpairlist,int *allpairlen)
 {
 int i,j,k,l,id1,id2,*assign;
 double dist,distsq,distcutoffsq;

 distcutoffsq = distcutoff*distcutoff;


 for(i=0;i<pointnum;i++)
  pairlen[i] = 0;

 for(i=0;i<featurenum;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   id1 = targetlist[i][j];
   for(k=j+1;k<targetlen[i];k++)
    {
    id2 = targetlist[i][k];
    distsq = getdistsq(coor[id1],coor[id2],DIM);
    if((distsq < distcutoffsq)&&(id2 > id1+3))
     {
     pairlist[id1][pairlen[id1]] = id2;
     pairlen[id1] ++;
     }
    }
   }
  }

 assign = intarray(pointnum);
 for(i=0;i<pointnum;i++)
  assign[i] = 1;

 getneighborlist_list_all_split_3(coor,pointnum,assign,allpairlist,allpairlen,distcutoff);

 free(assign);
 }
 





double getforce_feature(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **pairlist,int *pairlen,double **forcelist,double sphereradius)
 {
 int i,j,k,l,tempid;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF;

 contactcutoffsq = contactcutoff*contactcutoff;
 tempvector = doublearray(DIM);

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
//get nonbonded interaction
 for(i=0;i<pointnum;i++)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   distsq = getdistsq(coor[i],coor[tempid],DIM);
   if(distsq < contactcutoffsq)
    {
    dist = sqrt(distsq);
     
    for(k=0;k<DIM;k++)
     tempvector[k] = coor[tempid][k]-coor[i][k];
  
    force = 4*epsi*(-12*pow(sigma/dist,13)+6*pow(sigma/dist,7));

//    tempF = force;
//    if(tempF < 0)
//     tempF = -tempF;
//    if(tempF > tempFmax)
//     tempFmax = tempF;

    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] += force*tempvector[k]/dist;
     forcelist[tempid][k] -= force*tempvector[k]/dist;
     }
    }
   }
  }


//printf("vdw max %e\n",tempFmax);



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

//printf("Fmax %e\n",sqrt(Fmax));
//printf("\n");
 free(tempvector);
 return(sqrt(Fmax));
 }









double getforce_feature_all_noneighbor_pair_soft_norep_attraction(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen)
 {
 int i,j,k,l,tempid,id1,id2;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.2;
 tempvector = doublearray(DIM);
 factor = 0.0;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  tempid = pairlist[i][0];
  if((tempid > i+3)||(tempid < i-3))
   {
   distsq = getdistsq(coor[tempid],coor[i],DIM);
   distsq = sqrt(distsq);
   for(l=0;l<DIM;l++)
    tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
 
   force = epsi*(sigma/distsq-1);
//   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;


   for(l=0;l<DIM;l++)
    {
    forcelist[i][l] -= force*tempvector[l];
    forcelist[tempid][l] += force*tempvector[l];
    }
   }
  }
 }






 for(i=0;i<pointnum;i++)
  {
  for(j=i+4;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
   if(distsq < sigmasq)
    {
    distsq = sqrt(distsq);
    for(k=0;k<DIM;k++)
     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//    force = factor*epsi*(sigma/distsq-1);
    force = epsi*(sigma/distsq-1);
//    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   }
  }



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 



double getforce_feature_all_noneighbor_pair_soft_norep(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen)
 {
 int i,j,k,l,tempid,id1,id2;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*4.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  tempid = pairlist[i][0];
  if((tempid > i+1)||(tempid < i-1))
   {
   distsq = getdistsq(coor[tempid],coor[i],DIM);
   distsq = sqrt(distsq);
   for(l=0;l<DIM;l++)
    tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
 
   force = epsi*(sigma/distsq-1);
//   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;


   for(l=0;l<DIM;l++)
    {
    forcelist[i][l] -= force*tempvector[l];
    forcelist[tempid][l] += force*tempvector[l];
    }
   }
  }
 }






 for(i=0;i<pointnum;i++)
  {
  for(j=i+4;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    distsq = sqrt(distsq);
    for(k=0;k<DIM;k++)
     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
    force = factor*epsi*(sigma/distsq-1);
//    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   }
  }



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 



double getforce_feature_all_noneighbor_pair_soft(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen)
 {
 int i,j,k,l,tempid,id1,id2;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   if((tempid != i))
 //  if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    if(distsq > sigma)
     force = REPFAC*epsi*(sigma/distsq-1);
    else
     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }






 for(i=0;i<pointnum;i++)
  {
  for(j=i+2;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    distsq = sqrt(distsq);
    for(k=0;k<DIM;k++)
     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//    if(distsq > sigma)
//     force = factor*1*epsi*(sigma/distsq-1);
//    else
    if(distsq < sigma)
     force = REPFAC*epsi*(sigma/distsq-1);
//    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   }
  }



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 


double getforce_feature_all_noneighbor_pair_soft_sepattrep(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  avebondlen += dist;
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
//   if((tempid != i))
   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
//    force = REPFAC*(sigma/(distsq)-1)-epsi;
    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
//    if(distsq > sigma)
//     force = REPFAC*epsi*(sigma/distsq-1);
//    else
//     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }




testcount = 0;

 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    testcount ++;
//    distsq = sqrt(distsq);
    for(k=0;k<DIM;k++)
//     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
     tempvector[k] = (coor[j][k]-coor[i][k]);
//    if(distsq > sigma)
//     force = factor*1*epsi*(sigma/distsq-1);
//    else
    if(distsq < sigma)
    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
//     force = REPFAC*(sigma/(distsq)-1);
//     force = REPFAC*(sigma/(distsq)-1);

//     force = REPFAC*epsi*(sigma/distsq-1);
//    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   else
    {
    distsq = sqrt(distsq);
    j += (distsq-contactcutoff)/(avebondlen*1.2);
    }
   }
  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 



double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,int **contacttime)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  avebondlen += dist;
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
//    force = REPFAC*(sigma/(distsq)-1)-epsi;
//    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
//    dist6 = distsq*distsq*distsq;
    
    force = -epsi/sigma;
//    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
//    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
    force *= contacttime[i][j];
//    if(distsq > sigma)
//     force = REPFAC*epsi*(sigma/distsq-1);
//    else
//     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }



//for dark point
for(i=0;i<pointnum;i++)
 {
 if(randompairlen[i] != 0)
  {
  for(j=0;j<randompairlen[i];j++)
   {
   tempid = randompairlist[i][j];
//   if((tempid != i))
   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
//    force = REPFAC*(sigma/(distsq)-1)-epsi;
//    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
//    dist6 = distsq*distsq*distsq;
    force = -epsi/sigma;
//    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
//    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
//    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
//    if(distsq > sigma)
//     force = REPFAC*epsi*(sigma/distsq-1);
//    else
//     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }




testcount = 0;

 for(i=0;i<pointnum;i++)
  {
  for(j=i+1;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    testcount ++;
//    distsq = sqrt(distsq);
    for(k=0;k<DIM;k++)
//     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
     tempvector[k] = (coor[j][k]-coor[i][k]);
//    if(distsq > sigma)
//     force = factor*1*epsi*(sigma/distsq-1);
//    else
//    if(distsq < sigmasq)
    force = 0;
    if(distsq < sigmasq)
     {
     dist6 = distsq*distsq*distsq;
     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
     }
//    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
//     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
//     force = REPFAC*(sigma/(distsq)-1);
//     force = REPFAC*(sigma/(distsq)-1);

//     force = REPFAC*epsi*(sigma/distsq-1);
//    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   else
    {
    distsq = sqrt(distsq);
    j += (distsq-contactcutoff)/(avebondlen*1.2);
    }
   }
  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign,double *density,double *mass)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12,exfactor,exfactor6,exfactor12,shifta,tempdist,tempvalue;

 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 shifta = 3;
 exfactor = 0.2;

 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
   distsq = sqrt(distsq);
   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/distsq;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];

   if((tempid > i))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6)/distsq;
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }
    
    }
   }
  }
 }





//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  force /= mass[i];
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass_comp(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign,double *density,double *mass,int **compassign)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12,exfactor,exfactor6,exfactor12,shifta,tempdist,tempvalue;

 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 shifta = 3;
 exfactor = 0.2;

 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
   distsq = sqrt(distsq);
   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/distsq;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
//get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];

   if((tempid > i))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6)/distsq;
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }


//    if((compassign[i][1] != 1)&&(compassign[tempid][1] != 1)&&(distsq < 2*sigma))
//     {
//
//     for(l=0;l<DIM;l++)
//      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//
//     force = -10*epsi*(1/distsq-1/(2*sigma));
//
//     for(l=0;l<DIM;l++)
//      {
//      forcelist[i][l] -= force*tempvector[l];
//      forcelist[tempid][l] += force*tempvector[l];
//      }
//
//     }
    
    }
   }
  }
 }





//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  force /= mass[i];
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign,double *density)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12,exfactor,exfactor6,exfactor12,shifta,tempdist,tempvalue;

// WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 shifta = 3;
 exfactor = 0.2;
// exfactor = 1+shifta*pow(0.5,0.166666);
// exfactor6 = exfactor*exfactor*exfactor*exfactor*exfactor*exfactor;
// exfactor12 = exfactor6*exfactor6;

 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
   distsq = sqrt(distsq);
   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/distsq;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid > i))
//   if((tempid > i+3))
////   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
//    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
//    
//    force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   else if((tempid > i+1))
   if((tempid > i))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
//    if(((tempid==367)&&(i==370))||((tempid==370)&&(i==367)))
//     printf("id %d %d %lf %lf\n",i,tempid,distsq,WCAcutoff);
//    if(distsq < 0.3)
//     {
//     printf("id %d %d %lf %lf\n",i,tempid,distsq,WCAcutoff);
//     }
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }
    
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




//testcount = 0;
//for(i=0;i<pointnum;i++)
// {
// if(pairlen[i] != 0)
//  {
//  for(j=0;j<pairlen[i];j++)
//   {
//   tempid = pairlist[i][j];
////   if((tempid != i))
//
//
//   if((tempid > i))
////   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
//     for(l=0;l<DIM;l++)
//      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//   
//     tempvalue = exp(-exfactor*distsq);
//     force = -density[i]*density[tempid]*SCALE*epsi*exfactor*tempvalue/((1+tempvalue)*(1+tempvalue));
////     tempdist = distsq+sigma*shifta;
////     dist6 = tempdist*tempdist*tempdist*tempdist*tempdist*tempdist;
////     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
//     
// //    force = density[i]*density[tempid]*SCALE*4*epsi*(12*sigma12/(4096*dist6*dist6)-6*sigma6/(64*dist6));
////     force = density[i]*density[tempid]*SCALE*4*epsi*(12*exfactor6*exfactor6*sigma12/(dist6*dist6)-6*exfactor6*sigma6/(dist6));
//  
//  
//     for(l=0;l<DIM;l++)
//      {
//      forcelist[i][l] -= force*tempvector[l];
//      forcelist[tempid][l] += force*tempvector[l];
//      }
//    }
//
//
//
////   if((tempid > i))
////    {
////    distsq = getdistsq(coor[tempid],coor[i],DIM);
//////    distsq = sqrt(distsq);
////    for(l=0;l<DIM;l++)
////     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//// 
////    force = 0;
////    if(distsq < sigmasq)
////     {
////     dist6 = distsq*distsq*distsq;
////     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
////     }
//// 
////    for(l=0;l<DIM;l++)
////     {
////     forcelist[i][l] -= force*tempvector[l];
////     forcelist[tempid][l] += force*tempvector[l];
////     }
////    }
//
//
//
//   }
//  }
// }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12;

// WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
//  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid > i))
   if((tempid > i+3)||(tempid < i-3))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   else if((tempid > i+1)||(tempid < i-1))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
//    if(((tempid==367)&&(i==370))||((tempid==370)&&(i==367)))
//     printf("id %d %d %lf %lf\n",i,tempid,distsq,WCAcutoff);
//    if(distsq < 0.3)
//     {
//     printf("id %d %d %lf %lf\n",i,tempid,distsq,WCAcutoff);
//     }
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }
    
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




testcount = 0;
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
//   if((tempid != i))


   if((tempid > i+3))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = SCALE*4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }



//   if((tempid > i))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
// 
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }



   }
  }
 }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v2_density(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign,double *density)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12;

// WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
//  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid > i))
   if((tempid > i+3))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   else if(tempid > i+1)
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }
    
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




testcount = 0;
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
//   if((tempid != i))


   if((tempid > i))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = density[i]*density[tempid]*4*epsi*(12*sigma12/(4096*dist6*dist6)-6*sigma6/(64*dist6));
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }






   }
  }
 }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v2(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **allpairlist,int *allpairlen,int *chrassign)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6,sigma6,sigma12;

// WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoff = sigma*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma6 = sigma3*sigma3;
 sigma12 = sigma6*sigma6;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
//  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid > i))
   if((tempid > i+3))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   else if(tempid > i+1)
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    if((distsq<WCAcutoff))
     {
     for(l=0;l<DIM;l++)
      tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
   
     dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
     
     force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
  
  
     for(l=0;l<DIM;l++)
      {
      forcelist[i][l] -= force*tempvector[l];
      forcelist[tempid][l] += force*tempvector[l];
      }
     }
    
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




testcount = 0;
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
//   if((tempid != i))


   if((tempid > i+3))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
    dist6 = distsq*distsq*distsq*distsq*distsq*distsq;
    
    force = 4*epsi*(12*sigma12/(dist6*dist6)-6*sigma6/dist6);
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }






   }
  }
 }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 





 


double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,int **contacttime,int **allpairlist,int *allpairlen,int *chrassign)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  if(chrassign[i]==chrassign[i+1])
   {
   dist = getdistsq(coor[i],coor[i+1],DIM);
   dist = sqrt(dist);
   avebondlen += dist;
   for(j=0;j<DIM;j++)
    tempvector[j] = coor[i+1][j]-coor[i][j];
 
   force = 2*kb*(dist-bondlen);
 
 //  tempF = force;
 //  if(tempF < 0)
 //   tempF = -tempF;
 //  if(tempF > tempFmax)
 //   tempFmax = tempF;
 
 
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j] += force*tempvector[j]/dist;
    forcelist[i+1][j] -= force*tempvector[j]/dist;
    }
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
//    force = REPFAC*(sigma/(distsq)-1)-epsi;
//    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
//    dist6 = distsq*distsq*distsq;
    
    force = -0.01*epsi/sigma;
//    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
//    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
//    force *= contacttime[i][j];
//    if(distsq > sigma)
//     force = REPFAC*epsi*(sigma/distsq-1);
//    else
//     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




testcount = 0;
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid != i))
   if((tempid != i))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
 
    force = 0;
    if(distsq < sigmasq)
     {
     dist6 = distsq*distsq*distsq;
     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
     }
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 



double getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen,int **randompairlist,int *randompairlen,int **contacttime,int **allpairlist,int *allpairlen)
 {
 int i,j,k,l,tempid,id1,id2,testcount,temp,avebondlen;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq,sigma3,sigma5,sigma11,dist6;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*1.0;
 tempvector = doublearray(DIM);
 factor = 1.0;
 sigmasq = sigma*sigma;
 sigma3 = sigma*sigma*sigma;
 sigma5 = sigma*sigma*sigma*sigma*sigma;
 sigma11 = sigma5*sigma5*sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;
avebondlen = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  avebondlen += dist;
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);
avebondlen /= pointnum-1;


tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
//     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
  
//    force = REPFAC*(sigma/(distsq)-1)-epsi;
//    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
//    dist6 = distsq*distsq*distsq;
    
    force = -0.01*epsi/sigma;
//    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
//    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
//    force *= contacttime[i][j];
//    if(distsq > sigma)
//     force = REPFAC*epsi*(sigma/distsq-1);
//    else
//     force = REPFAC*epsi*(sigma/distsq-1);
     
 //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }



////for dark point
//for(i=0;i<pointnum;i++)
// {
// if(randompairlen[i] != 0)
//  {
//  for(j=0;j<randompairlen[i];j++)
//   {
//   tempid = randompairlist[i][j];
////   if((tempid != i))
//   if((tempid > i+3)||(tempid < i-3))
//    {
//    distsq = getdistsq(coor[tempid],coor[i],DIM);
////    distsq = sqrt(distsq);
//    for(l=0;l<DIM;l++)
//     tempvector[l] = (coor[tempid][l]-coor[i][l]);
////     tempvector[l] = (coor[tempid][l]-coor[i][l])/distsq;
//  
////    force = REPFAC*(sigma/(distsq)-1)-epsi;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-RANDOMFAC*epsi/sigma;
////    dist6 = distsq*distsq*distsq;
//    force = -epsi/sigma;
////    force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6))-epsi/sigma;
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq))-epsi/sigma;
////    force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq)-epsi/sigma;
////    if(distsq > sigma)
////     force = REPFAC*epsi*(sigma/distsq-1);
////    else
////     force = REPFAC*epsi*(sigma/distsq-1);
//     
// //   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
// 
// 
//    for(l=0;l<DIM;l++)
//     {
//     forcelist[i][l] -= force*tempvector[l];
//     forcelist[tempid][l] += force*tempvector[l];
//     }
//    }
//   }
//  }
// }




testcount = 0;
for(i=0;i<pointnum;i++)
 {
 if(allpairlen[i] != 0)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
//   if((tempid != i))
   if((tempid != i))
    {
    distsq = getdistsq(coor[tempid],coor[i],DIM);
//    distsq = sqrt(distsq);
    for(l=0;l<DIM;l++)
     tempvector[l] = (coor[tempid][l]-coor[i][l]);
 
    force = 0;
    if(distsq < sigmasq)
     {
     dist6 = distsq*distsq*distsq;
     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
     }
 
    for(l=0;l<DIM;l++)
     {
     forcelist[i][l] -= force*tempvector[l];
     forcelist[tempid][l] += force*tempvector[l];
     }
    }
   }
  }
 }

// for(i=0;i<pointnum;i++)
//  {
//  for(j=i+1;j<pointnum;j++)
//   {
//   distsq = getdistsq(coor[i],coor[j],DIM);
//   if(distsq < contactcutoffsq)
//    {
//    testcount ++;
////    distsq = sqrt(distsq);
//    for(k=0;k<DIM;k++)
////     tempvector[k] = (coor[j][k]-coor[i][k])/distsq;
//     tempvector[k] = (coor[j][k]-coor[i][k]);
////    if(distsq > sigma)
////     force = factor*1*epsi*(sigma/distsq-1);
////    else
////    if(distsq < sigmasq)
//    force = 0;
//    if(distsq < sigmasq)
//     {
//     dist6 = distsq*distsq*distsq;
//     force = REPFAC*(sigma11/(dist6*dist6)-sigma5/(dist6));
//     }
////    force = REPFAC*(sigma5/(distsq*distsq*distsq)-sigma3/(distsq*distsq));
////     force = REPFAC*(sigma3/(distsq*distsq)-sigma/distsq);
////     force = REPFAC*(sigma/(distsq)-1);
////     force = REPFAC*(sigma/(distsq)-1);
//
////     force = REPFAC*epsi*(sigma/distsq-1);
////    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//    for(k=0;k<DIM;k++)
//     {
//     forcelist[i][k] -= force*tempvector[k];
//     forcelist[j][k] += force*tempvector[k];
//     }
//    }
//   else
//    {
//    distsq = sqrt(distsq);
//    j += (distsq-contactcutoff)/(avebondlen*1.2);
//    }
//   }
//  }


//printf("testcount %d\n",testcount);

 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

 free(tempvector);
 return(sqrt(Fmax));
 }




 




double getforce_feature_all_noneighbor_pair(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius,int **pairlist,int *pairlen)
 {
 int i,j,k,l,tempid,id1,id2;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff*0.8;
 tempvector = doublearray(DIM);
 factor = 0.05;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
////get nonbonded interaction
for(i=0;i<pointnum;i++)
 {
 if(pairlen[i] != 0)
  {
  tempid = pairlist[i][0];
  if((tempid > i+3)||(tempid < i-3))
   {
   distsq = getdistsq(coor[tempid],coor[i],DIM);
   for(l=0;l<DIM;l++)
    tempvector[l] = coor[tempid][l]-coor[i][l];
 
   force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;


   for(l=0;l<DIM;l++)
    {
    forcelist[i][l] -= force*tempvector[l];
    forcelist[tempid][l] += force*tempvector[l];
    }
   }
  }
 }

// for(i=0;i<featurenum;i++)
//  {
//  for(j=0;j<targetlen[i];j++)
//   {
//   for(k=j+1;k<targetlen[i];k++)
//    {
//    id1 = targetlist[i][j];
//    id2 = targetlist[i][k];
//    if((id1 > id2+3)||(id1 < id2-3))
//     {
//     distsq = getdistsq(coor[id1],coor[id2],DIM);
//     if(distsq < 4*sigmasq)
//      {
//       
//      for(l=0;l<DIM;l++)
//       tempvector[l] = coor[id2][l]-coor[id1][l];
//    
//      force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
//  
//  
//      for(l=0;l<DIM;l++)
//       {
//       forcelist[id1][l] -= force*tempvector[l];
//       forcelist[id2][l] += force*tempvector[l];
//       }
//      }
//     }
//    }
//   }
//  }


 for(i=0;i<pointnum;i++)
  {
  for(j=i+4;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    for(k=0;k<DIM;k++)
     tempvector[k] = coor[j][k]-coor[i][k];
    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   }
  }



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

//printf("Fmax %e %e %e\n",sqrt(Fmax),sigmasq,contactcutoffsq);
//printf("Fmax %e\n",sqrt(Fmax));
//printf("\n");
 free(tempvector);
 return(sqrt(Fmax));
 }




 




double getforce_feature_all_noneighbor(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **targetlist,int featurenum,int *targetlen,double **forcelist,double sphereradius)
 {
 int i,j,k,l,tempid,id1,id2;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff;
 tempvector = doublearray(DIM);
 factor = 0.05;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
//get nonbonded interaction
 for(i=0;i<featurenum;i++)
  {
  for(j=0;j<targetlen[i];j++)
   {
   for(k=j+1;k<targetlen[i];k++)
    {
    id1 = targetlist[i][j];
    id2 = targetlist[i][k];
    if((id1 > id2+3)||(id1 < id2-3))
     {
     distsq = getdistsq(coor[id1],coor[id2],DIM);
     if(distsq < 4*sigmasq)
      {
       
      for(l=0;l<DIM;l++)
       tempvector[l] = coor[id2][l]-coor[id1][l];
    
      force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
  
  
      for(l=0;l<DIM;l++)
       {
       forcelist[id1][l] -= force*tempvector[l];
       forcelist[id2][l] += force*tempvector[l];
       }
      }
     }
    }
   }
  }


 for(i=0;i<pointnum;i++)
  {
  for(j=i+4;j<pointnum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < contactcutoffsq)
    {
    for(k=0;k<DIM;k++)
     tempvector[k] = coor[j][k]-coor[i][k];
    force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
    for(k=0;k<DIM;k++)
     {
     forcelist[i][k] -= force*tempvector[k];
     forcelist[j][k] += force*tempvector[k];
     }
    }
   }
  }



 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

//printf("Fmax %e %e %e\n",sqrt(Fmax),sigmasq,contactcutoffsq);
//printf("Fmax %e\n",sqrt(Fmax));
//printf("\n");
 free(tempvector);
 return(sqrt(Fmax));
 }




 




double getforce_feature_all(double **coor,int pointnum,double bondlen,double kb,double contactcutoff,double sigma,double epsi,int **pairlist,int *pairlen,double **forcelist,double sphereradius,int **allpairlist,int *allpairlen)
 {
 int i,j,k,l,tempid;
 double distsq,dist,*tempvector,coef,force,tempcutoff,Fmax,contactcutoffsq,tempFmax,tempF,WCAcutoff,WCAcutoffsq,factor,sigmasq;

 WCAcutoff = contactcutoff*pow(2,1/6.0);
 WCAcutoffsq = WCAcutoff*WCAcutoff;
 contactcutoffsq = contactcutoff*contactcutoff;
 tempvector = doublearray(DIM);
 factor = 0.05;
 sigmasq = sigma*sigma;

 for(i=0;i<pointnum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

tempFmax = 0;
 //getsphere restriction
 for(i=0;i<pointnum;i++)
  {
  distsq = getnormsq(coor[i],DIM);
  if(distsq > sphereradius*sphereradius)
   {
//   tempF = sqrt(distsq/sphereradius);
//   if(tempF > tempFmax)
//    tempFmax = tempF;

   for(j=0;j<DIM;j++)
    forcelist[i][j] = -coor[i][j]/sphereradius;
   }
  }

//printf("sphere max %e\n",tempFmax);

tempFmax = 0;

 //get bonded interaction
 for(i=0;i<pointnum-1;i++)
  {
  dist = getdistsq(coor[i],coor[i+1],DIM);
  dist = sqrt(dist);
  for(j=0;j<DIM;j++)
   tempvector[j] = coor[i+1][j]-coor[i][j];

  force = 2*kb*(dist-bondlen);

//  tempF = force;
//  if(tempF < 0)
//   tempF = -tempF;
//  if(tempF > tempFmax)
//   tempFmax = tempF;


  for(j=0;j<DIM;j++)
   {
   forcelist[i][j] += force*tempvector[j]/dist;
   forcelist[i+1][j] -= force*tempvector[j]/dist;
   }
  }

//printf("bond max %e\n",tempFmax);



tempFmax = 0;
//get nonbonded interaction
 for(i=0;i<pointnum;i++)
  {
  for(j=0;j<pairlen[i];j++)
   {
   tempid = pairlist[i][j];
   if((tempid > i+3)||(tempid < i-3))
    {
    distsq = getdistsq(coor[i],coor[tempid],DIM);
//    if(distsq < contactcutoffsq)
    if(distsq < 4*contactcutoffsq)
     {
//     dist = sqrt(distsq);
      
     for(k=0;k<DIM;k++)
      tempvector[k] = coor[tempid][k]-coor[i][k];
   
     force = 4*epsi*(pow(sigmasq/distsq,6.5)-pow(sigmasq/distsq,3.5))/sigma;
 
//     tempF = force;
//     if(tempF < 0)
//      tempF = -tempF;
//     if(tempF > tempFmax)
//      tempFmax = tempF;
 
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] -= force*tempvector[k];
      forcelist[tempid][k] += force*tempvector[k];
      }
     }
    }
   }
  }



//printf("vdw max %e\n",tempFmax);


tempFmax = 0;
//get nonbonded interaction for all pair
 for(i=0;i<pointnum;i++)
  {
  for(j=0;j<allpairlen[i];j++)
   {
   tempid = allpairlist[i][j];
   if((i>tempid+3)||(i<tempid-3))
    {
    distsq = getdistsq(coor[i],coor[tempid],DIM);
    if(distsq < contactcutoffsq)
     {
//     dist = sqrt(distsq);
      
     for(k=0;k<DIM;k++)
      tempvector[k] = coor[tempid][k]-coor[i][k];
  
     force = factor*4*epsi*(pow(contactcutoffsq/distsq,6.5)-pow(contactcutoffsq/distsq,3.5))/contactcutoff;
//     force = factor*4*epsi*(-12*pow(sigma/dist,13)+6*pow(sigma/dist,7));

//     tempF = force;
//     if(tempF < 0)
//      tempF = -tempF;
//     if(tempF > tempFmax)
//      tempFmax = tempF;

     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] -= force*tempvector[k];
      forcelist[tempid][k] += force*tempvector[k];
      }
     }
    }

   }
  }

//printf("all vdw max %e\n",tempFmax);





 Fmax = 0;
 for(i=0;i<pointnum;i++)
  {
  force = getnormsq(forcelist[i],DIM);
  if(force > Fmax)
   Fmax = force;
  }

//printf("Fmax %e %e %e\n",sqrt(Fmax),sigmasq,contactcutoffsq);
printf("Fmax %e\n",sqrt(Fmax));
printf("\n");
 free(tempvector);
 return(sqrt(Fmax));
 }





void shrinkcoor_v2(double **coor,int statenum,double shrinkfactor)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] *= shrinkfactor;
 }






double  getmaxdist(double **coor, int pointnum)
 {
 int i,j,k;
 double maxdist,dist;

 maxdist = 0;
 for(i=0;i<pointnum;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   dist += coor[i][j]*coor[i][j];
  if(dist > maxdist)
   maxdist = dist;
  }

 maxdist = sqrt(maxdist);
 return(maxdist);

 }






void main()
{
char inputfilename[N],assignfilename[N],outputfilename[N],densityfilename[N],massfilename[N],compfilename[N];
int i,j,k,l,m,n,pointnum,featurenum,*assign,**featureassign,stepnum,interval,seed,assignmax,**targetlist,*targetlen,**pairlist,*pairlen,pairinterval,shrinkmax,shrinkinterval,shrinkendbool,**randompairlist,*randompairlen,sepprob,**contacttime,**allpairlist,*allpairlen,*chrassign,tempmax,maxlen,*tempassign,nextstepnum,outputnum,**compassign;
double **coor,**velvector,**velvector2,**forcelist,radius,kT,bondlen,kb,beadradius,Fmax,dt,dtp,t,sphereradius,chrgyr,pointmass,sigma,epsi,distcutoff,factor,gyr_target,*density,densitybase,center[DIM],maxdist,*mass,cutoffsphereradius;
FILE *inputfile,*assignfile,*outputfile,*densityfile,*massfile,*compfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the assignment for chromatin point:\n");
scanf("%s",assignfilename);

//printf("Input the number of chr bead in the system:\n");
//scanf("%d",&chrlen);
//
printf("Input the filename for density:\n");
scanf("%s",densityfilename);

printf("Input the base for density normalization:\n");
scanf("%lf",&densitybase);

printf("Input filename for mass:\n");
scanf("%s",massfilename);

printf("Input filename for compartment:\n");
scanf("%s",compfilename);

printf("Input the value for bond lenth:\n");
scanf("%lf",&bondlen);

//printf("Input the kb for bond vibration:\n");
//scanf("%lf",&kb);
//
printf("Input the step number for iteration:\n");
scanf("%d",&stepnum);

printf("Input the step number for iteration of the second round:\n");
scanf("%d",&nextstepnum);

printf("Input the interval for saving:\n");
scanf("%d",&interval);

printf("Input the value for random number generator:\n");
scanf("%d",&seed);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


srand(seed);
kT=1;
sepprob = 0.1;
beadradius = bondlen/2;
distcutoff = beadradius*4;
sigma = beadradius*2.0;
//epsi = 1.0*kT;	//WANG CHEN XI 
epsi = 0.5*kT;	//WANG CHEN XI 
//nearcutoff = 100;
//nearcutoff = bondlen;
//contactprob = 0.1;
//sepprob = 0.1;
kb = 100*kT/(bondlen*bondlen);
dtp = bondlen/10;
pointmass = 1;
pairinterval = distcutoff*0.5/dtp;
if(pairinterval < 1)
 pairinterval = 1;
 
//shrinkmax = 1000;
shrinkmax = 100000;
shrinkinterval = 100;
maxlen = 200;	//set by mldaet

outputnum = 100;
//interval = 20000;
//interval = (stepnum+nextstepnum)/(2*outputnum);


pointnum = getlinenum(inputfilename);

assign = intarray(pointnum);
tempassign = intarray(pointnum);
coor = doublematrix(pointnum,DIM);
chrassign = intarray(pointnum);
velvector = doublematrix(pointnum,DIM);
velvector2 = doublematrix(pointnum,DIM);
forcelist = doublematrix(pointnum,DIM);
density = doublearray(pointnum);
mass = doublearray(pointnum);
compassign = intmatrix(pointnum,2);

for(i=0;i<DIM;i++)
 center[i] = 0;

inputfile = openfile(inputfilename,'r');
massfile = openfile(massfilename,'r');
for(i=0;i<pointnum;i++)
 {
 fscanf(massfile,"%lf",&mass[i]);
 tempassign[i] = 1;
 for(j=0;j<DIM;j++)
  {
  fscanf(inputfile,"%lf",&coor[i][j]);
  center[j] += coor[i][j];
  }
 }
fclose(inputfile);
fclose(massfile);

for(i=0;i<DIM;i++)
 {
 center[i] /= pointnum;
 for(j=0;j<pointnum;j++)
  coor[j][i] -= center[i];
 }


compfile = openfile(compfilename,'r');
for(i=0;i<pointnum;i++)
 for(j=0;j<2;j++)
  fscanf(compfile,"%d",&compassign[i][j]);
fclose(compfile);


assignmax = -1;
assignfile = openfile(assignfilename,'r');
for(i=0;i<pointnum;i++)
 {
 chrassign[i] = compassign[i][0];
// fscanf(assignfile,"%d",&chrassign[i]);
 fscanf(assignfile,"%d",&assign[i]);
 if(assignmax < assign[i])
  assignmax = assign[i];
 }
fclose(assignfile);

densityfile = openfile(densityfilename,'r');
for(i=0;i<pointnum;i++)
 {
 fscanf(densityfile,"%lf",&density[i]);
 density[i] /= densitybase;
 }
fclose(densityfile);

featurenum = (int)(log((double)assignmax)/log(2))+1;

featureassign = intmatrix(pointnum,featurenum);

for(i=0;i<pointnum;i++)
 {
 getfeatureassign(featureassign[i],featurenum,assign[i]);
// getfeatureassign_randomtype(featureassign[i],featurenum+1,assign[i]);
 }

targetlist = intmatrix(featurenum,pointnum);
targetlen = intarray(featurenum);
//pairlist = intpointarray(pointnum);
pairlist = intmatrix(pointnum,maxlen);
//pairlist = intmatrix(pointnum,pointnum*(featurenum));
pairlen = intarray(pointnum);
//randompairlist = intmatrix(pointnum,pointnum);
//randompairlen = intarray(pointnum);
//contacttime = intmatrix(pointnum,pointnum);
allpairlist = intmatrix(pointnum,maxlen);
//allpairlist = intmatrix(pointnum,maxlen);
allpairlen = intarray(pointnum);

gettargetlist(featureassign,pointnum,featurenum,targetlist,targetlen);

sphereradius = beadradius*pow(pointnum/0.35,0.3333);
//cutoffsphereradius = beadradius*pow(pointnum,0.3333)*3.0;
cutoffsphereradius = sphereradius;
gyr_target = sphereradius*0.77;
//gyr_target = 80;

shrinkendbool = 0;
maxlen = 100000;
//getpairlist_mon(coor,targetlist,targetlen,featurenum,pointnum,pairlist,pairlen,2*sigma);
//getpairlist_mon(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,1.5*sigma);

outputfile = openfile(outputfilename,'w');
for(i=0;i<stepnum+nextstepnum;i++)
 {
// genvel(velvector,pointnum,DIM,kT);
 genvel_mass(velvector,pointnum,DIM,kT,mass);


// printf("Step %d %d\n",i,shrinkendbool);
 if((i%shrinkinterval == 0)&&(shrinkendbool == 0)&&(i<shrinkmax))
  {
//  if(i==0)
   shrinkcoor_v2(coor,pointnum,0.99);
   maxdist = getmaxdist(coor,pointnum);
   if(maxdist < cutoffsphereradius)
    shrinkendbool = 1;
//   chrgyr = getgyr(coor,pointnum);
//   if(chrgyr < gyr_target)
//    shrinkendbool = 1;
  }
// if(i%interval == 0)
// if((i>stepnum+nextstepnum-outputnum*interval)&&(i%interval==0))
 if((i%interval==0))
  {
   chrgyr = getgyr(coor,pointnum);
  printf("%d %lf %lf %lf %e %d\n",i,t,chrgyr,sphereradius,Fmax,shrinkendbool);
//  fprintf(outputfile,"Snapshot %d %lf %lf\n",i,t,chrgyr);
//  fprintf(outputfile,"Snapshot %d %lf %lf\n",i,t,Fmax);
//  for(j=0;j<pointnum;j++)
//   fprintf(outputfile,"%lf %lf %lf\n",coor[j][0],coor[j][1],coor[j][2]);
  }

 if(i%10==0)
  {
//  getpairlist_mon(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,1.2*sigma);
//  getpairlist_random(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,randompairlist,randompairlen,1.2*sigma);
//  getpairlist_random_allagg(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,randompairlist,randompairlen,1.2*sigma,assign,featureassign);
//  getpairlist_random_allagg_prob(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,randompairlist,randompairlen,1.2*sigma,assign,featureassign);
//  getpairlist_random_allagg_prob_time(coor,targetlist,targetlen,featurenum+1,pointnum,pairlist,pairlen,randompairlist,randompairlen,1.4*sigma,assign,featureassign,contacttime);
//  getpairlist_random_allagg_prob_time_alllist(coor,targetlist,targetlen,featurenum,pointnum,pairlist,pairlen,randompairlist,randompairlen,2.0*sigma,assign,featureassign,contacttime,allpairlist,allpairlen);
//  getpairlist_random_allagg_prob_time_alllist_v2(coor,targetlist,targetlen,featurenum,pointnum,pairlist,pairlen,2.0*sigma,assign,featureassign,allpairlist,allpairlen,tempassign);
  getpairlist_random_allagg_prob_time_alllist_v2(coor,targetlist,targetlen,featurenum,pointnum,pairlist,pairlen,4.0*beadradius,assign,featureassign,allpairlist,allpairlen,tempassign);
//  getpairlist_random_allagg_prob_time_alllist_v3(coor,targetlist,targetlen,featurenum,pointnum,pairlist,pairlen,2.0*sigma,sphereradius,assign,featureassign,allpairlist,allpairlen,tempassign);
// printf("check dist %d %d %lf\n",367,370,getdistsq(coor[367],coor[370],DIM));
  tempmax = allpairlen[0];
//  for(j=0;j<allpairlen[367];j++)
//   printf("contact %d %d\n",j,allpairlist[367][j]);
  for(j=0;j<pointnum;j++)
   if(tempmax < allpairlen[j])
    tempmax = allpairlen[j];
//  if(maxlen <= tempmax)
//   {
//   printf("Error for the maxlen setting\n");
//   exit(0);
//   }
//  else
   chrgyr = getgyr(coor,pointnum);
   if(i%1000==0)
    {
    printf("maxlen %d\n",tempmax);
    printf("%d %lf %lf %lf %e %d\n",i,t,chrgyr,sphereradius,Fmax,shrinkendbool);
    }
//  printf("%d %lf %lf %lf %lf %e\n",i,t,chrgyr,sphereradius,beadradius,Fmax);
  }

// Fmax=getforce_feature_all_noneighbor(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius);
// Fmax=getforce_feature_all_noneighbor_pair(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen);
// Fmax=getforce_feature_all_noneighbor_pair_soft(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen);
// Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum+1,targetlen,forcelist,sphereradius,pairlist,pairlen,randompairlist,randompairlen,contacttime);
// Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,randompairlist,randompairlen,contacttime,allpairlist,allpairlen);
// Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,randompairlist,randompairlen,contacttime,allpairlist,allpairlen,chrassign);
 if(i < stepnum)
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v2(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign);
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density);
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density,mass);
  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass_comp(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density,mass,compassign);
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v2_density(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density);
 else
  {
  epsi = 0.5*kT;

//  if(((i-stepnum)%shrinkinterval == 0)&&((i-stepnum)<300))
//   {
//   shrinkcoor_v2(coor,pointnum,0.99);
//   }
 

// printf("check dist %d %d %d %lf\n",i,367,370,getdistsq(coor[367],coor[370],DIM));

  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass_comp(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density,mass,compassign);
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density_mass(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density,mass);
//  Fmax=getforce_feature_all_noneighbor_pair_soft_sepattrep_darkatt_allpair_chr_v3_density(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen,allpairlist,allpairlen,chrassign,density);
  }
  
// Fmax=getforce_feature_all_noneighbor_pair_soft_norep_attraction(coor,pointnum,bondlen,kb,beadradius*2,sigma,epsi,targetlist,featurenum,targetlen,forcelist,sphereradius,pairlist,pairlen);

 if((i == 0)&&(Fmax < 1))
  Fmax = 1;
 dt = sqrt(2*dtp/Fmax);

// printf("check %d %lf %lf\n",i,t,dt);
// if(i%10==0)
//  printf("%d %lf %lf %lf %e %d\n",i,t,chrgyr,sphereradius,Fmax,shrinkendbool);
//  printf("Step %d Fmax %e %lf\n",i,Fmax,dt);
 t += dt;


// getfinalvel(velvector,pointnum,forcelist,dt,velvector2);
 getfinalvel_mass(velvector,pointnum,forcelist,dt,velvector2,mass);
// printf("vel %lf %lf %lf %lf\n\n",velvector[0][0],forcelist[0][0],dt,velvector2[0][0]);
 updatecoor(coor,velvector,velvector2,pointnum,dt);
 
 }


//fprintf(outputfile,"Snapshot %d %lf %lf\n",i,t,chrgyr);
fprintf(outputfile,"Snapshot %d %lf %lf\n",i,t,Fmax);
for(j=0;j<pointnum;j++)
 fprintf(outputfile,"%lf %lf %lf\n",coor[j][0],coor[j][1],coor[j][2]);

fclose(outputfile);




for(i=0;i<featurenum;i++)
 free(targetlist[i]);
free(targetlist);

for(i=0;i<pointnum;i++)
 {
 free(pairlist[i]);
 free(coor[i]);
 free(velvector[i]);
 free(velvector2[i]);
 free(forcelist[i]);
 free(featureassign[i]);
// free(randompairlist[i]);
// free(contacttime[i]);
 free(allpairlist[i]);
 free(compassign[i]);
 }

free(allpairlist);
//free(contacttime);
free(targetlen);
free(pairlen);
free(assign);
//free(randompairlen);
free(allpairlen);

free(pairlist);
free(coor);
free(velvector);
free(velvector2);
free(forcelist);
free(featureassign);
//free(randompairlist);
free(chrassign);
//free(tempassign);
free(compassign);



}










