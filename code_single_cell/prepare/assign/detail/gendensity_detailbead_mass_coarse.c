#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
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
char pointlenfilename[N],outputpointlenfilename[N],assignfilename[N],compassignfilename[N],outputcompassignfilename[N];
int i,j,id,k,ratio,statenum_old,statenum_new,*pointlen_old,*pointlen_new,*assign,totalpointnum,*compartassign_old,*compartassign_new;
FILE *pointlenfile,*outputpointlenfile,*assignfile,*compassignfile,*outputcompassignfile;


printf("Input the filename for original pointlen:\n");
scanf("%s",pointlenfilename);

printf("Input the filename for compartment assignment:\n");
scanf("%s",compassignfilename);

printf("Input the ratio for coarse-graining:\n");
scanf("%d",&ratio);

printf("Input the filename for output pointlen:\n");
scanf("%s",outputpointlenfilename);

printf("Input the filename for output assignment:\n");
scanf("%s",assignfilename);

printf("Input the filename for output compartment assignment:\n");
scanf("%s",outputcompassignfilename);

statenum_old = getlinenum(pointlenfilename);
statenum_new = (statenum_old-1)/ratio+1;

if(statenum_old != getlinenum(compassignfilename))
 {
 printf("ERROR in the input file:\n");
 exit(0);
 }


pointlen_old = intarray(statenum_old);
pointlen_new = intarray(statenum_new);
compartassign_old = intarray(statenum_old);
compartassign_new = intarray(statenum_new);

pointlenfile = openfile(pointlenfilename,'r');
compassignfile = openfile(compassignfilename,'r');
for(i=0;i<statenum_old;i++)
 {
 fscanf(pointlenfile,"%d",&pointlen_old[i]);
 fscanf(compassignfile,"%d",&compartassign_old[i]);
 }
fclose(pointlenfile);
fclose(compassignfile);

for(i=0;i<statenum_new;i++)
 {
 pointlen_new[i] = 0;
 compartassign_new[i] = 0;
 }

for(i=0;i<statenum_old;i++)
 {
 id = i/ratio;
 pointlen_new[id] += pointlen_old[i];
 if(compartassign_old[i] == 1)
  compartassign_new[id] ++;
 }

totalpointnum = 0;
for(i=0;i<statenum_new;i++)
 {
 pointlen_new[i] /= ratio;
 if(pointlen_new[i] == 0)
  pointlen_new[i] = 1;
 totalpointnum += pointlen_new[i];
 if(compartassign_new[i] > ratio/2.0001)
  compartassign_new[i] = 1;
 else
  compartassign_new[i] = 0;
 }

assign = intarray(totalpointnum);

for(i=0;i<totalpointnum;i++)
 assign[i] = -1;

id=0;
for(i=0;i<statenum_new;i++)
 {
 for(j=0;j<pointlen_new[i];j++)
  {
  assign[id] = i;
  id ++;
  }
 }


outputpointlenfile = openfile(outputpointlenfilename,'w');
outputcompassignfile = openfile(outputcompassignfilename,'w');
for(i=0;i<statenum_new;i++)
 {
 fprintf(outputpointlenfile,"%d\n",pointlen_new[i]);
 fprintf(outputcompassignfile,"%d\n",compartassign_new[i]);
 }
fclose(outputpointlenfile);
fclose(outputcompassignfile);

assignfile = openfile(assignfilename,'w');
for(i=0;i<totalpointnum;i++)
 fprintf(assignfile,"%d\n",assign[i]);
fclose(assignfile);

free(pointlen_old);
free(pointlen_new);
free(assign);
free(compartassign_old);
free(compartassign_new);


}


