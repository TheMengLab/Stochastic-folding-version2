#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3


void main()
{
char inputfilename[N],assignfilename[N],outputfilename[N];
int i,j,k,linenum,pointnum,*assign,*pointlen;
double **coor,**mergecoor;
FILE *inputfile,*assignfile,*outputfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

linenum = getlinenum(inputfilename);

assign = intarray(linenum);
coor = doublematrix(linenum,DIM);

inputfile = openfile(inputfilename,'r');
assignfile = openfile(assignfilename,'r');
for(i=0;i<linenum;i++)
 {
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
 fscanf(assignfile,"%d",&assign[i]);
 }
fclose(inputfile);
fclose(assignfile);

pointnum = assign[linenum-1]+1;

pointlen = intarray(pointnum);
mergecoor = doublematrix(pointnum,DIM);

for(i=0;i<pointnum;i++)
 for(j=0;j<DIM;j++)
  mergecoor[i][j] = 0;

for(i=0;i<pointnum;i++)
 pointlen[i] = 0;

for(i=0;i<linenum;i++)
 {
 pointlen[assign[i]] ++;
 for(j=0;j<DIM;j++)
  mergecoor[assign[i]][j] += coor[i][j];
 }

for(i=0;i<pointnum;i++)
 for(j=0;j<DIM;j++)
  mergecoor[i][j] /= pointlen[i];

outputfile = openfile(outputfilename,'w');
for(i=0;i<pointnum;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",mergecoor[i][j]);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);

for(i=0;i<linenum;i++)
 free(coor[i]);
free(coor);

for(i=0;i<pointnum;i++)
 {
 free(mergecoor[i]);
 }
free(mergecoor);
free(assign);
free(pointlen);




}



