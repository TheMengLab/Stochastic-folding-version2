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






double getscore(double **scorematrix,int statenum,int bid,int eid)
 {
 int i,j,k,count1,count2;
 double value1,value2,score;

 count1 = 0;
 count2 = 0;
 value1 = 0;
 value2 = 0;

 for(i=bid;i<eid;i++)
  {
  for(j=i;j<i+eid-bid;j++)
   {
   if((i>=0)&&(i<statenum)&&(j>=0)&&(j<statenum))
    {
    if(j-bid > 2*(i-bid))
     {
     value1 += scorematrix[i][j];
     count1 ++;
     }
    else
     {
     value2 += scorematrix[i][j];
     count2 ++;
     }
    }
   }
  }

 score = 0;

 if((count1!=0)&&(count2!=0))
  score = value1/count1-value2/count2;

 return(score);


 }





double getscore_assign(double **scorematrix,int statenum,int bid,int eid,int *assign)
 {
 int i,j,k,count1,count2;
 double value1,value2,score;

 count1 = 0;
 count2 = 0;
 value1 = 0;
 value2 = 0;

 for(i=bid;i<eid;i++)
  {
  for(j=i;j<i+eid-bid;j++)
   {
   if((i>=0)&&(i<statenum)&&(j>=0)&&(j<statenum))
    {
    if(j-bid > 2*(i-bid))
     {
     if(2*i-j > 0)
      {
      if((assign[i] == 1)&&(assign[j] == 1)&&(assign[2*i-j] == 1))
       {
       value1 += scorematrix[i][j];
       count1 ++;
       }
      }
     }
    else
     {
     if(2*i-j > 0)
      {
      if((assign[i] == 1)&&(assign[j] == 1)&&(assign[2*i-j] == 1))
       {
       value2 += scorematrix[i][j];
       count2 ++;
       }
      }
     
     }
    }
   }
  }

 score = -1;

 if((count1!=0)&&(count2!=0))
  score = value1/count1-value2/count2;

 return(score);


 }







void main()
{
char inputfilename[N],outputfilename[N],scorefilename[N],assignfilename[N];
int i,j,k,statenum,linenum,*assign,bid,eid,lag,minlag,count,assignlen,*compassign;
double **matrix,**scorematrix,*scorelist,tempscore;
FILE *inputfile,*outputfile,*scorefile,*assignfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the dimension of matrix:\n");
scanf("%d",&statenum);

printf("Input the filename for compartassign(1 col):\n");
scanf("%s",assignfilename);

printf("Input filename for arrowhead output:\n");
scanf("%s",outputfilename);

printf("Input the filename for score list:\n");
scanf("%s",scorefilename);


assignlen = getlinenum(assignfilename);

matrix = doublematrix(statenum,statenum);
scorematrix = doublematrix(statenum,statenum);
assign = intarray(statenum);
scorelist = doublearray(statenum);
compassign = intarray(assignlen);

lag = 30;
minlag = 10;


assignfile = openfile(assignfilename,'r');
for(i=0;i<assignlen;i++)
 fscanf(assignfile,"%d",&compassign[i]);
fclose(assignfile);


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
 count = 0;
 for(j=0;(j<statenum);j++)
  {
  if(matrix[i][j] > 0.01)
   count ++;
//   assign[i] = 1;
  }
// if(count > statenum*0.01)
  assign[i] = 1;
 scorelist[i] = 0;
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


for(i=0;i<statenum;i++)
 {
 if(assign[i] == 1)
  {
  bid = i;
  eid = i+lag;

  if(eid > statenum-1)
   eid = statenum-1;

  scorelist[i] = getscore_assign(scorematrix,statenum,bid,eid,assign);
  }
 }


for(i=0;i<statenum;i++)
 {
 for(j=minlag;j<lag;j++)
  {
  if(assign[i] == 1)
   {
   bid = i;
   eid = i+j;
 
   if(eid > statenum-1)
    eid = statenum-1;
 
   tempscore = getscore_assign(scorematrix,statenum,bid,eid,assign);
   if(tempscore > scorelist[i])
    scorelist[i] = tempscore;
   }
  }
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



scorefile = openfile(scorefilename,'w');
for(i=0;i<statenum;i++)
 {
//  fprintf(outputfile,"%lf\n",scorelist[i]);
 if(i>=assignlen)
  fprintf(outputfile,"%lf\n",scorelist[i]);
 else if(compassign[i] != -1)
  fprintf(outputfile,"%lf\n",scorelist[i]);
 else
  fprintf(outputfile,"-1\n");

 }
fclose(scorefile);




for(i=0;i<statenum;i++)
 {
 free(matrix[i]);
 free(scorematrix[i]);
 }
free(matrix);
free(scorematrix);
free(assign);
free(scorelist);
free(compassign);


}








