#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3




void main()
{
char boundaryfilename[N],Tboundaryfilename[N],outputfilename[N];
int i,j,k,boundarynum,Tboundarynum,*boundary;
double temp,*Tboundary,dist,mindist;
FILE *boundaryfile,*Tboundaryfile,*outputfile;


printf("Input filename for boundary:\n");
scanf("%s",boundaryfilename);

printf("Input filename for TAD boundary:\n");
scanf("%s",Tboundaryfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


boundarynum = getlinenum(boundaryfilename);
Tboundarynum = getlinenum(Tboundaryfilename);

boundary = intarray(boundarynum);
Tboundary = doublearray(Tboundarynum);


boundaryfile = openfile(boundaryfilename,'r');
for(i=0;i<boundarynum;i++)
 {
 fscanf(boundaryfile,"%d",&boundary[i]);
 fscanf(boundaryfile,"%lf",&temp);
 }
fclose(boundaryfile);


Tboundaryfile = openfile(Tboundaryfilename,'r');
for(i=0;i<Tboundarynum;i++)
 {
 fscanf(Tboundaryfile,"%lf",&Tboundary[i]);
 fscanf(Tboundaryfile,"%lf",&temp);
 }
fclose(Tboundaryfile);


outputfile = openfile(outputfilename,'w');
for(i=0;i<boundarynum;i++)
 {
 mindist = 1000000;
 for(j=0;j<Tboundarynum;j++)
  {
  dist = boundary[i]-Tboundary[j];
  if(dist < 0)
   dist = -dist;
  if(dist < mindist)
   mindist = dist;
  }
 fprintf(outputfile,"%d %lf\n",boundary[i],mindist);
 }
fclose(outputfile);

free(boundary);
free(Tboundary);

}





