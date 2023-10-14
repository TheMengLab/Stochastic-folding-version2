#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#define N 200
#define DIM 3
#define RATIO 5



void main()
{
char matrix1filename[N],matrix2filename[N],outputfilename[N];
int i,j,k,statenum,data1,data2,coarsenum,id1,id2;
double **matrix,**coarsematrix,value;
FILE *matrix1file,*matrix2file,*outputfile;

printf("Input filename for the matrix:\n");
scanf("%s",matrix1filename);

printf("How many states in this system:\n");
scanf("%d",&statenum);

printf("Input filename foru output matrix:\n");
scanf("%s",outputfilename);


coarsenum = (int)(statenum-1)/RATIO+1;

//matrix = doublematrix(statenum,statenum);
coarsematrix = doublematrix(coarsenum,coarsenum);

for(i=0;i<coarsenum;i++)
 for(j=0;j<coarsenum;j++)
  coarsematrix[i][j] = 0;

matrix1file = openfile(matrix1filename,'r');
for(i=0;i<statenum;i++)
 {
 id1 = i/RATIO;
 for(j=0;j<statenum;j++)
  {
  id2 = j/RATIO;
  fscanf(matrix1file,"%lf",&value);
  coarsematrix[id1][id2] += value;
  }
 }
fclose(matrix1file);

outputfile = openfile(outputfilename,'w');

for(i=0;i<coarsenum;i++)
 {
 for(j=0;j<coarsenum;j++)
  fprintf(outputfile,"%lf ",coarsematrix[i][j]);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);


for(i=0;i<coarsenum;i++)
 free(coarsematrix[i]);
free(coarsematrix);


}



