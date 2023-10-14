#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/matrix/inverse/matrixinverse.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
#define N 500
#define DIM 3
#define FACTOR 1
#define MAXD 0.005

void getcoef(double *ylist,double **Xmatrix,int ynum,int xnum,double *coef)
 {
 int i,j,k;
 double **tempmatrix,**tempmatrix_inverse;


 tempmatrix = doublematrix(xnum,xnum);
 tempmatrix_inverse = doublematrix(xnum,xnum);

 for(i=0;i<xnum;i++)
  for(j=0;j<xnum;j++)
   {
   tempmatrix[i][j] = 0;
   for(k=0;k<ynum;k++)
    tempmatrix[i][j] += Xmatrix[k][i]*Xmatrix[k][j];
   }
 
 matrixinverse_lapack(tempmatrix, xnum, tempmatrix_inverse);


 for(i=0;i<xnum;i++)
  {
  coef[i] = 0;
  for(j=0;j<ynum;j++)
   for(k=0;k<xnum;k++)
    coef[i] += ylist[j]*Xmatrix[j][k]*tempmatrix_inverse[i][k];
  }

 for(i=0;i<xnum;i++)
  {
  free(tempmatrix[i]);
  free(tempmatrix_inverse[i]);
  }
 free(tempmatrix);
 free(tempmatrix_inverse);
 }






double getdistsq(double *coor1,double *coor2,int dim)
 {
 int i,j;
 double distsq;


 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(distsq);
 }







void getassign(int *pointlen,int statenum,int *assign,int linenum,int *startid)
 {
 int i,j,k,id;

 for(i=0;i<linenum;i++)
  assign[i] = -1;
 id = 0;
 for(i=0;i<statenum;i++)
  {
//  printf("%d %d %d\n",i,pointlen[i],id);
//  if(pointlen[i] > 20)
//   printf("%d %d\n",i,pointlen[i]);
  for(j=0;j<pointlen[i];j++)
   {
   if(j==0)
    startid[i] = id;

   assign[id] = i;
   id ++;
   if(id > linenum)
    {
    printf("ERROR when do the assignment:\n");
    exit(0);
    }
   }
  }
 }





void gettempcoor(double ***coor,int filenum,int linenum,double ***tempcoor,int *assign,int *pointlen,int *startid,int statenum)
 {
 int i,j,k,l,len;

 for(i=0;i<filenum;i++)
  {
  for(j=0;j<statenum;j++)
   {
   len = pointlen[j];
   for(k=0;k<DIM;k++)
    {
    tempcoor[i][j][k] = 0;
    for(l=0;l<len;l++)
     tempcoor[i][j][k] += coor[i][startid[j]+l][k];
    tempcoor[i][j][k] /= len;
    }
   }
  }

 }






double getmediandist(double ***tempcoor,int filenum,int statenum,int *targetpair,double *distlist)
 {
 int i,j,k,l,id1,id2;
 double dist,tempvalue;

 id1 = targetpair[0];
 id2 = targetpair[1];

 for(i=0;i<filenum;i++)
  {
  distlist[i] = getdistsq(tempcoor[i][id1],tempcoor[i][id2],DIM);
  }

 for(i=0;i<=filenum/2;i++)
  {
  for(j=i+1;j<filenum;j++)
   {
   if(distlist[i] < distlist[j])
    {
    tempvalue = distlist[i];
    distlist[i] = distlist[j];
    distlist[j] = tempvalue;
    }
   }

  }
 dist = log(distlist[filenum/2])/2;

 return(dist);
 }





double getavedist(double ***tempcoor,int filenum,int statenum,int *targetpair)
 {
 int i,j,k,l,id1,id2;
 double dist,tempvalue;

 id1 = targetpair[0];
 id2 = targetpair[1];

 dist = 0;
 for(i=0;i<filenum;i++)
  {
  tempvalue = sqrt(getdistsq(tempcoor[i][id1],tempcoor[i][id2],DIM));
//  distlist[i] = sqrt(getdistsq(tempcoor[i][id1],tempcoor[i][id2],DIM));
  dist += tempvalue;
  }
 dist /= filenum;

// for(i=0;i<=filenum/2;i++)
//  {
//  for(j=i+1;j<filenum;j++)
//   {
//   if(distlist[i] < distlist[j])
//    {
//    tempvalue = distlist[i];
//    distlist[i] = distlist[j];
//    distlist[j] = tempvalue;
//    }
//   }
//
//  }
 dist = log(dist);

 return(dist);
 }





double getroundvalue(double data)
 {
 int intpart;
 double decpart,targetdata;

 intpart = (int)data;
 if(data < 0)
  intpart --;
 decpart = data - intpart;

 if(rand()/((double)RAND_MAX) > decpart)
  targetdata = intpart;
 else
  targetdata = intpart +1;
 return(targetdata);

 }



void getshiftmaxminlist(double *shift,int statenum,int winsize,int *maxidlist,int *minidlist)
 {
 int i,j,k,startid,endid;

 for(i=0;i<statenum;i++)
  {
  maxidlist[i] = 1;
  minidlist[i] = 1;
  }


 for(i=0;i<statenum;i++)
  {
  startid = i-winsize;
  endid = i+winsize;

  if(startid < 0)
   startid = 0;
  if(endid > statenum-1)
   endid = statenum-1;

  for(j=startid;j<=endid;j++)
   {
   if(shift[j] > shift[i])
    maxidlist[i] = 0;
   if(shift[j] < shift[i])
    minidlist[i] = 0;
   }

  }


 }





void getfitcontactvalue(double *pointlen,int statenum,int **contactidlist,double *fitcontactvalue,int contactnum)
 {
 int i,j,k,id1,id2;
 double dist,contactvalue;

 for(i=0;i<contactnum;i++)
  {
  id1 = contactidlist[i][0];
  id2 = contactidlist[i][1];
  dist = (pointlen[id1]+pointlen[id2])/2;

  for(j=id1+1;j<id2;j++)
   dist += pointlen[j];

  contactvalue = -log(dist);

  fitcontactvalue[i] = contactvalue;
  }


 }







double getcontactdiff(double *contactvalue,double *fitcontactvalue,int contactnum,double ave1,double ave2)
 {
 int i,j,k;
 double diff,shift;


 diff = 0;
 shift = ave2-ave1;

 for(i=0;i<contactnum;i++)
  {
  diff += (fitcontactvalue[i]-contactvalue[i]-shift)*(fitcontactvalue[i]-contactvalue[i]-shift);
  }

 return(sqrt(diff/contactnum));
 }






double getavevalue(double *contactvalue,int contactnum)
 {
 int i,j;
 double ave;

 ave = 0;
 for(i=0;i<contactnum;i++)
  ave += contactvalue[i];
 ave /= contactnum;

 return(ave);

 }







void getpointlenshiftvalue(double *pointlen,int statenum,int **contactidlist,double *contactvalue,double *fitcontactvalue,double ave1,double ave2,double *shift,int contactnum)
 {
 int i,j,k,id1,id2,len;
 double targetdist,targetcontact,dist,maxdiff,diff;

 for(i=0;i<statenum;i++)
  shift[i] = 0;


 for(i=0;i<contactnum;i++)
  {
  id1 = contactidlist[i][0];
  id2 = contactidlist[i][1];

  targetcontact = contactvalue[i]-ave1+ave2;
  targetdist = 1/exp(targetcontact);

  dist = 0;
  dist = (pointlen[id1]+pointlen[id2])/2;

  for(j=id1+1;j<id2;j++)
   dist += pointlen[j];

  diff = targetdist-dist;
  if(diff > MAXD)
   diff = MAXD;
  if(diff < -MAXD)
   diff = -MAXD;
  shift[id1] += 0.5*diff/(id2-id1);
  shift[id2] += 0.5*diff/(id2-id1);
  for(j=id1+1;j<id2;j++)
   shift[j] += diff/(id2-id1);
  }

 for(i=0;i<statenum;i++)
  {
  pointlen[i] += shift[i];
  if(pointlen[i] < 1)
   pointlen[i] = 1;
  if(pointlen[i] > 20)
   pointlen[i] = 20;
  }


 }





void getproblist(double lambda,double *problist,int len)
 {
 int i,j,k;

 problist[0] = exp(-lambda);
 for(i=1;i<len;i++)
  {
  problist[i] = problist[i-1]*lambda/i;
  }
 }



int getpossionlen(double *problist,int len,double tempprob)
 {
 int i,j,id,endbool;

 id = -1;
 endbool = 0;
 for(i=0;i<len;i++)
  {
  if(tempprob < problist[i])
   {
   id = i;
   endbool = 1;
   }
  else
   tempprob -= problist[i];
  }

 if(id == -1)
  id = len;

 return(id);
 }







double getseqdist(double *len,int id1,int id2)
 {
 int i,j;
 double dist;

 dist = (len[id1]+len[id2])/2;
 for(i=id1+1;i<id2;i++)
  dist += len[i];

 return(dist);

 }







double getfitcontactvalue_matrix(double **len,int repnum,int id1,int id2)
 {
 int i,j,k;
 double dist,fitvalue,contactvalue;

 contactvalue=0;
 for(i=0;i<repnum;i++)
  {
  dist = getseqdist(len[i],id1,id2);
  contactvalue += 1/dist;
  }
 
 contactvalue /= repnum;

 return(contactvalue);
 }







void main()
{
char inputfilename[N],outputfilename[N],pointlenfilename[N];
int i,j,k,l,m,statenum,repnum,seed,nonzeronum,range,minrange,minrange2;
double *data,*tempdata,basevalue,shiftvalue,**len,**matrix,database,*problist,tempprob,*rowsum;
FILE *inputfile,*outputfile,*pointlenfile;

printf("Input the filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the seed for random number generator:\n");
scanf("%d",&seed);

printf("Input the value for open base:\n");
scanf("%lf",&basevalue);

printf("Input the value for open shift:\n");
scanf("%lf",&shiftvalue);

printf("Input the value for replica num:\n");
scanf("%d",&repnum);

printf("Input the filename for output matrix:\n");
scanf("%s",outputfilename);

printf("Input the filename for output pointlen:\n");
scanf("%s",pointlenfilename);


//repnum = 20;
range = 100;
minrange = 1;
minrange2 = 5;

statenum = getlinenum(inputfilename);
srand(seed);

data = doublearray(statenum);
tempdata = doublearray(statenum);
len = doublematrix(repnum,statenum);
//matrix = doublematrix(statenum,statenum);
problist = doublearray(100);
//rowsum = doublearray(statenum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 fscanf(inputfile,"%lf",&data[i]);
fclose(inputfile);


nonzeronum = 0;

for(i=0;i<statenum;i++)
 {
 if(data[i] > 0.00001)
  {
  tempdata[nonzeronum] = data[i];
  nonzeronum ++;
  }
// rowsum[i] =0;
 }

database = getmedianvalue(tempdata,nonzeronum);

//for(i=0;i<statenum;i++)
// data[i] /= database;

//get random length
for(j=0;j<statenum;j++)
 {
// getproblist(data[j]/basevalue,problist,100);
 getproblist(data[j]/database,problist,100);
 for(i=0;i<repnum;i++)
  {
  tempprob = rand()/((double)(RAND_MAX));
  len[i][j] = basevalue+shiftvalue*getpossionlen(problist,100,tempprob);
  }
 }

//for(i=0;i<statenum;i++)
// for(j=0;j<statenum;j++)
//  {
//  matrix[i][j] = 0;
//  }


//for(i=0;i<statenum;i++)
// {
// for(j=i+1;j<statenum;j++)
//  {
//  if((j-i<range)&&(j-i>minrange))
//   {
//   matrix[i][j] = getfitcontactvalue_matrix(len,repnum,i,j);
//   matrix[j][i] = matrix[i][j];
//   }
//  }
// }

//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  rowsum[i] += matrix[i][j];
// }


//outputfile = openfile(outputfilename,'w');
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  if((j-i > minrange2)||(j-i < -minrange2))
//   fprintf(outputfile,"%e ",matrix[i][j]/sqrt(rowsum[i]*rowsum[j]));
//  else
//   fprintf(outputfile,"0 ");
//   
// fprintf(outputfile,"\n");
// }
//fclose(outputfile);



pointlenfile = openfile(pointlenfilename,'w');
for(i=0;i<repnum;i++)
 {
 for(j=0;j<statenum;j++)
  fprintf(pointlenfile,"%d ",(int)(len[i][j]+0.5));
 fprintf(pointlenfile,"\n");
 }
fclose(pointlenfile);



free(data);
free(tempdata);
free(problist);
for(i=0;i<repnum;i++)
 free(len[i]);
free(len);
//for(i=0;i<statenum;i++)
// free(matrix[i]);
//free(matrix);
//free(rowsum);

}








