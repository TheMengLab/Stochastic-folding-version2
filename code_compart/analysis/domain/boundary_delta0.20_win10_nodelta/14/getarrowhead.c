#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
//#define R_SPHERE_SQ 90000








//void gen3Drandomvector(double *vector)
// {
// int i,j;
// double theta,phi;
//
// theta = rand()/(double)RAND_MAX*2-1;
// phi = rand()/(double)RAND_MAX*2*PI;
//
// theta = acos(theta);
//
// vector[0] = sin(theta)*cos(phi);
// vector[1] = sin(theta)*sin(phi);
// vector[2] = cos(theta);
// }




double getdistsq(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;
 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(distsq);
 }






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





void main()
{
char inputfilename[N],outputfilename[N];
int i,j,k,statenum,linenum,*assign;
double **matrix,**scorematrix;
FILE *inputfile,*outputfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the dimension of matrix:\n");
scanf("%d",&statenum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


matrix = doublematrix(statenum,statenum);
scorematrix = doublematrix(statenum,statenum);
assign = intarray(statenum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  {
  fscanf(inputfile,"%lf",&matrix[i][j]);
  scorematrix[i][j] = 0;
  }
fclose(inputfile);



for(i=0;i<statenum;i++)
 {
 assign[i] = 0;
 for(j=0;(j<statenum)&&(assign[i]==0);j++)
  {
  if(matrix[i][j] > 0.01)
   assign[i] = 1;
  }
 }

for(i=0;i<statenum;i++)
 for(j=i+1;j<statenum;j++)
  {
  if((2*i-j > 0)&&(assign[i]==1)&&(assign[j]==1))
   {
   if(matrix[i][j]+matrix[i][2*i-j] > 0.01)
    scorematrix[i][j] = (matrix[i][j]-matrix[i][2*i-j])/(matrix[i][j]+matrix[i][2*i-j]);
   }
  scorematrix[j][i] = scorematrix[i][j];
  }


outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  {
  fprintf(outputfile,"%lf ",scorematrix[i][j]);
  }
 fprintf(outputfile,"\n");
 }
fclose(outputfile);


for(i=0;i<statenum;i++)
 {
 free(matrix[i]);
 free(scorematrix[i]);
 }
free(matrix);
free(scorematrix);
free(assign);


}








