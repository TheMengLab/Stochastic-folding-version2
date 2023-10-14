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
int i,j,id,k,ratio,statenum_old,statenum_new,*pointlen_old,*pointlen_new,*assign,totalpointnum,**compartassign_old,**compartassign_new,*chrlen_old,*chrlen_new,chrnum,base1,base2,id1,id2;
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

if(statenum_old != getlinenum(compassignfilename))
 {
 printf("ERROR in the input file:\n");
 exit(0);
 }

pointlen_old = intarray(statenum_old);
compartassign_old = intmatrix(statenum_old,2);

chrnum = 0;
pointlenfile = openfile(pointlenfilename,'r');
compassignfile = openfile(compassignfilename,'r');
for(i=0;i<statenum_old;i++)
 {
 fscanf(pointlenfile,"%d",&pointlen_old[i]);
 for(j=0;j<2;j++)
  fscanf(compassignfile,"%d",&compartassign_old[i][j]);
 if(chrnum<compartassign_old[i][0])
  chrnum = compartassign_old[i][0];
 }
fclose(pointlenfile);
fclose(compassignfile);

chrnum ++;

chrlen_old = intarray(chrnum);
chrlen_new = intarray(chrnum);

for(i=0;i<chrnum;i++)
 {
 chrlen_old[i] = 0;
 chrlen_new[i] = 0;
 }

for(i=0;i<statenum_old;i++)
 {
 chrlen_old[compartassign_old[i][0]] ++;
 }

statenum_new = 0;
for(i=0;i<chrnum;i++)
 {
 chrlen_new[i] = (chrlen_old[i]-1)/ratio+1;
 statenum_new += chrlen_new[i];
 }


pointlen_new = intarray(statenum_new);
compartassign_new = intmatrix(statenum_new,2);


for(i=0;i<statenum_new;i++)
 {
 pointlen_new[i] = 0;
 compartassign_new[i][0] = 0;
 compartassign_new[i][1] = 0;
 }

base1 = 0;
base2 = 0;
for(i=0;i<chrnum;i++)
 {
 for(j=0;j<chrlen_old[i];j++)
  {
  id1 = j+base1;
  id2 = j/ratio+base2;
  if(id2 >= statenum_new)
   {
   printf("\n\nERROR IN CODE\n\n");
   exit(0);
   }
  pointlen_new[id2] += pointlen_old[id1];

  compartassign_new[id2][0] = i;
  if(compartassign_old[id1][1] == 1)
   compartassign_new[id2][1] ++;
  }
 base1 += chrlen_old[i];
 base2 += chrlen_new[i];
 }


totalpointnum = 0;
for(i=0;i<statenum_new;i++)
 {
 pointlen_new[i] /= ratio;
 if(pointlen_new[i] == 0)
  pointlen_new[i] = 1;
 totalpointnum += pointlen_new[i];

 if(compartassign_new[i][1] > ratio/1.999)
  compartassign_new[i][1] = 1;
 else if(compartassign_new[i][1] > ratio/2.0001)
  {
  if(rand()%10000 > 5000)
   compartassign_new[i][1] = 1;
  else
   compartassign_new[i][1] = 0;
  }
 else
  compartassign_new[i][1] = 0;
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
 fprintf(outputcompassignfile,"%d %d\n",compartassign_new[i][0],compartassign_new[i][1]);
 }
fclose(outputpointlenfile);
fclose(outputcompassignfile);

assignfile = openfile(assignfilename,'w');
for(i=0;i<totalpointnum;i++)
 fprintf(assignfile,"%d\n",assign[i]);
fclose(assignfile);









for(i=0;i<statenum_old;i++)
 free(compartassign_old[i]);
free(compartassign_old);

for(i=0;i<statenum_new;i++)
 free(compartassign_new[i]);
free(compartassign_new);

free(chrlen_old);
free(chrlen_new);




free(pointlen_old);
free(pointlen_new);
free(assign);


}


