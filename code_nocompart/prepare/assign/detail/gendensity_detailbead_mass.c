#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define PI 3.14159265358979323846264338327950288419716939937510
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
//#define R_SPHERE_SQ 90000


int getassign(int *dataassign,int linenum,int id,int ratio)
 {
 int i,j,assign;

 assign = 0;

 for(i=id*ratio;(i<id*ratio+ratio)&&(i<linenum);i++)
  if(dataassign[i]==1)
   assign = 1;

 return(assign);
 }




void main()
{
char inputfilename[N],outputfilename[N],pointlenfilename[N],massfilename[N];
int i,j,k,l,pointnum,*pointlen,*outputassign,totalpointnum,maxlen;
double *density,interval;
FILE *inputfile,*outputfile,*pointlenfile,*massfile;

printf("Input filename for density data:\n");
scanf("%s",inputfilename);

printf("Input filename for interval:\n");
scanf("%lf",&interval);

printf("Input the maxlen for each bin:\n");
scanf("%d",&maxlen);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

printf("Input filename for point len:\n");
scanf("%s",pointlenfilename);

printf("Input filename for mass:\n");
scanf("%s",massfilename);

pointnum = getlinenum(inputfilename);

pointlen = intarray(pointnum);
density = doublearray(pointnum);

totalpointnum = 0;
inputfile = openfile(inputfilename,'r');
for(i=0;i<pointnum;i++)
 {
 fscanf(inputfile,"%lf",&density[i]);
 pointlen[i] = (int)(density[i]/interval)+1;
 if(pointlen[i] > maxlen)
  pointlen[i] = maxlen;
 totalpointnum += pointlen[i];
 }
fclose(inputfile);


outputassign = intarray(totalpointnum);
j=0;
for(i=0;i<pointnum;i++)
 {
 for(k=0;k<pointlen[i];k++)
  {
  outputassign[j] = i;
  j ++;
//  if(pointlen[i] == 1)
//   printf("0\n");
//  else
//   printf("1\n");
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<totalpointnum;i++)
 {
 fprintf(outputfile,"%d\n",outputassign[i]);
 }
fclose(outputfile);

pointlenfile = openfile(pointlenfilename,'w');
for(i=0;i<pointnum;i++)
 fprintf(pointlenfile,"%d\n",pointlen[i]);
fclose(pointlenfile);

massfile = openfile(massfilename,'w');
for(i=0;i<pointnum;i++)
 {
 for(j=0;j<pointlen[i];j++)
  fprintf(massfile,"%lf\n",1.0/pointlen[i]);
 }
fclose(massfile);


free(outputassign);
free(pointlen);
free(density);



}


