#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getrmsd.c"
#define N 200
#define DIM 3




double getdistsq(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;

 distsq = 0;
 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(distsq);
 }





void main()
{
char inputlistfilename[N],inputfilename[N],outputfilename[N],targetfilename[N];
int i,j,k,id1,id2,filenum,atomnum,*targetid,*targetassign,targetnum,toplinenum;
double **matrix,***coor,distcutoff,distsq,distcutoffsq,*sum;
FILE *inputlistfile,*inputfile,*outputfile,*targetfile;

printf("Input the filename for filelist of original coordinate:\n");
scanf("%s",inputlistfilename);

printf("Input the number of atoms in this system:\n");
scanf("%d",&atomnum);

printf("Input filename for target atomlist(start from 1):\n");
scanf("%s",targetfilename);

printf("Input the cutoff for contact dist(with radius of 10):\n");
scanf("%lf",&distcutoff);

printf("Input the value for top linenum:\n");
scanf("%d",&toplinenum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


filenum = getlinenum(inputlistfilename);
targetnum = getlinenum(targetfilename);

targetid = intarray(targetnum);
targetassign = intarray(targetnum);
coor = doublematrixarray(filenum,atomnum,DIM);
sum = doublearray(targetnum);
matrix = doublematrix(targetnum,toplinenum);


inputlistfile = openfile(inputlistfilename,'r');
for(i=0;i<filenum;i++)
 {
 fscanf(inputlistfile,"%s",inputfilename);
 inputfile = openfile(inputfilename,'r');
 for(j=0;j<atomnum;j++)
  for(k=0;k<DIM;k++)
   fscanf(inputfile,"%lf",&coor[i][j][k]);
 fclose(inputfile);
 }
fclose(inputlistfile);

targetfile = openfile(targetfilename,'r');
for(i=0;i<targetnum;i++)
 {
 fscanf(targetfile,"%d",&targetid[i]);
 fscanf(targetfile,"%d",&targetassign[i]);
 targetid[i] --;
 }
fclose(targetfile);


for(i=0;i<targetnum;i++)
 {
 for(j=0;j<toplinenum;j++)
  matrix[i][j] = 0;
 sum[i] = 0;
 }

distcutoffsq = distcutoff*distcutoff;
printf("%lf %lf\n",distcutoff,distcutoffsq);
for(i=0;i<targetnum;i++)
 {
 id1 = targetid[i];


 for(j=1;j<toplinenum;j++)
  {
  if(i+j < targetnum)
   {
   id2 = targetid[i+j];
   if((targetassign[i] != -1)&&(targetassign[i+j] != -1))
    {
    for(k=0;k<filenum;k++)
     {
     distsq = getdistsq(coor[k][id1],coor[k][id2],DIM);
     if(distsq < distcutoffsq)
      matrix[i][j] ++;
     }
//    matrix[j][i] = matrix[i][j];
    }
   }


//  id2 = targetid[j];
//  if((targetassign[i] != -1)&&(targetassign[j] != -1))
//   {
//   for(k=0;k<filenum;k++)
//    {
//    distsq = getdistsq(coor[k][id1],coor[k][id2],DIM);
//    if(distsq < distcutoffsq)
//     matrix[i][j] ++;
//    }
//   matrix[j][i] = matrix[i][j];
//   }
  }


// for(j=i+1;j<targetnum;j++)
//  {
//  id2 = targetid[j];
//  if((targetassign[i] != -1)&&(targetassign[j] != -1))
//   {
//   for(k=0;k<filenum;k++)
//    {
//    distsq = getdistsq(coor[k][id1],coor[k][id2],DIM);
//    if(distsq < distcutoffsq)
//     matrix[i][j] ++;
//    }
//   matrix[j][i] = matrix[i][j];
//   }
//  }

 }


for(i=0;i<targetnum;i++)
 {
 id1 = targetid[i];

 for(j=i+1;j<targetnum;j++)
  {
  id2 = targetid[j];
  if((targetassign[i] != -1)&&(targetassign[j] != -1))
   {
   for(k=0;k<filenum;k++)
    {
    distsq = getdistsq(coor[k][id1],coor[k][id2],DIM);
    if(distsq < distcutoffsq)
     {
     sum[i] ++;
     sum[j] ++;
     }
    }
   }
  }

 }



printf("test\n");

outputfile = openfile(outputfilename,'w');
for(i=0;i<targetnum;i++)
 {
 for(j=0;j<toplinenum;j++)
  {
  if(i+j<targetnum)
   {
   if(sum[i]*sum[i+j] > 0.001)
    fprintf(outputfile,"%e ",matrix[i][j]/sqrt(sum[i]*sum[i+j])); 
   else
    fprintf(outputfile,"0 ");
   }
  else
   fprintf(outputfile,"0 ");
  }
// for(j=0;j<targetnum;j++)
//  fprintf(outputfile,"%e ",matrix[i][j]/sqrt(sum[i]*sum[j]));
 fprintf(outputfile,"\n");
 }
fclose(outputfile);





for(i=0;i<filenum;i++)
 {
 for(j=0;j<atomnum;j++)
  free(coor[i][j]);
 free(coor[i]);
 }
free(coor);

for(i=0;i<targetnum;i++)
 free(matrix[i]);
free(matrix);

free(targetid);
free(targetassign);
free(sum);


}




