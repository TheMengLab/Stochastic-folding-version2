#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3




int getendid(int *assign,int totallen,int targetid,int startid)
 {
 int i,j,k;

 int endid,endbool;

 endid = startid;

 endbool = 0;
 for(;endbool == 0;)
  {
  if(endid >= totallen)
   endbool = 1;
  else if(assign[endid] != targetid)
   {
   endbool = 1;
   }
  else
   endid ++;
  }

 endid --;

 return(endid);

 }





int getseglistid(double *seglist,int totallen,double refvalue)
 {
 int i,j,k,startid,endid,tempid,endbool;

 startid = 0;
 endid = totallen;

 for(;endid != startid+1;)
  {
  tempid = (startid+endid)/2;

  if(seglist[tempid] > refvalue+0.00001)
   {
   endid = tempid;
   }
  else if(seglist[tempid] < refvalue-0.00001)
   startid = tempid;
  else
   {
   startid = tempid;
   endid = startid+1;
   }
  }

 return(startid);

 }






void main()
{
char inputfilename[N],assign1filename[N],assign2filename[N],outputfilename[N],chrassign1filename[N],chrassign2filename[N];
int i,j,k,statenum,insertnum,chrnum,*assign1,*chrlen1,totallen1,*assign2,*chrlen2,ratio,totallen2,tempid1,tempid2,*chrassign1,*chrassign2,pointnum1,pointnum2,startid,endid,count,*chrbid1,*chrbid2;
double **coor,**coor2,*seg1list,*seg2list,tempdist,range,expratio,shiftvalue;
FILE *inputfile,*assign1file,*assign2file,*outputfile,*chrassign1file,*chrassign2file;

printf("Input the filename for the conformation:\n");
scanf("%s",inputfilename);

printf("Input filename for first assignment:\n");
scanf("%s",assign1filename);

printf("Input the filename for first chr assignment:\n");
scanf("%s",chrassign1filename);

//printf("Input the ratio for insert:\n");
//scanf("%d",&ratio);
//
printf("Input filename for second assignment:\n");
scanf("%s",assign2filename);

printf("Input the filename for second chr assignment:\n");
scanf("%s",chrassign2filename);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);


totallen1 = getlinenum(assign1filename);
totallen2 = getlinenum(assign2filename);
statenum = getlinenum(inputfilename);

pointnum1 = getlinenum(chrassign1filename);
pointnum2 = getlinenum(chrassign2filename);

seg1list = doublearray(totallen1);
seg2list = doublearray(totallen2);


chrnum = 1;

insertnum = (int)(pointnum2/((double)pointnum1)+0.5);
printf("%d\n",insertnum);
expratio = pow(insertnum,1/3.0);


if(totallen1 != statenum)
 {
 printf("ERROR for the input\n");
 exit(0);
 }

assign1 = intarray(totallen1);
assign2 = intarray(totallen2);
chrassign1 = intarray(pointnum1);
chrassign2 = intarray(pointnum2);

coor = doublematrix(statenum,DIM);
coor2 = doublematrix(totallen2,DIM);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
fclose(inputfile);


//printf("%lf %lf\n",coor[100][0],coor[500][0]);


assign1file = openfile(assign1filename,'r');
for(i=0;i<totallen1;i++)
 {
// fscanf(assign1file,"%d",&j);
 fscanf(assign1file,"%d",&assign1[i]);
// fscanf(assign1file,"%d",&j);
// chrlen1[assign1[i]] ++;
 }
fclose(assign1file);

assign2file = openfile(assign2filename,'r');
for(i=0;i<totallen2;i++)
 {
// fscanf(assign2file,"%d",&j);
 fscanf(assign2file,"%d",&assign2[i]);
// fscanf(assign2file,"%d",&j);
// chrlen2[assign2[i]] ++;
 }
fclose(assign2file);



chrassign1file = openfile(chrassign1filename,'r');
for(i=0;i<pointnum1;i++)
 {
 fscanf(chrassign1file,"%d",&chrassign1[i]);
 if(chrnum < chrassign1[i])
  chrnum = chrassign1[i];
 fscanf(chrassign1file,"%d",&j);
 }
fclose(chrassign1file);


chrassign2file = openfile(chrassign2filename,'r');
for(i=0;i<pointnum2;i++)
 {
 fscanf(chrassign2file,"%d",&chrassign2[i]);
 fscanf(chrassign2file,"%d",&j);
 }
fclose(chrassign2file);

chrnum ++;
chrlen1 = intarray(chrnum);
chrlen2 = intarray(chrnum);

chrbid1 = intarray(chrnum);
chrbid2 = intarray(chrnum);

for(i=0;i<chrnum;i++)
 {
 chrbid1[i] = 0;
 chrbid2[i] = 0;
 chrlen1[i] = 0;
 chrlen2[i] = 0;
 }


for(i=1;i<totallen1;i++)
 {
 if(chrassign1[assign1[i]] != chrassign1[assign1[i-1]])
  {
  tempid1 = chrassign1[assign1[i]];
  chrbid1[tempid1] = i;
  }
 }

for(i=1;i<totallen2;i++)
 {
 if(chrassign2[assign2[i]] != chrassign2[assign2[i-1]])
  {
  tempid1 = chrassign2[assign2[i]];
  chrbid2[tempid1] = i;
  }
 }


for(i=0;i<pointnum1;i++)
 {
 chrlen1[chrassign1[i]] ++;
 }

for(i=0;i<pointnum2;i++)
 {
 chrlen2[chrassign2[i]] ++;
 }






startid = 0;
//calculate seg list
for(i=0;i<pointnum1;i++)
 {
 if(assign1[startid] != i)
  {
  printf("Error in the assignment1:\n");
  printf("%d %d %d %d\n",startid,assign1[startid],i,pointnum1);
  exit(0);
  }
 else
  {
  endid = getendid(assign1,totallen1,i,startid);
  count = endid-startid+1;
  for(j=startid;j<=endid;j++)
   {
   seg1list[j] = i+(j-startid)*1.0/(count);
   seg1list[j] *= insertnum;
   }
  startid = endid+1;
  }
 }


startid = 0;
//calculate seg list
for(i=0;i<pointnum2;i++)
 {
 if(assign2[startid] != i)
  {
  printf("Error in the assignment2:\n");
  exit(0);
  }
 else
  {
  endid = getendid(assign2,totallen2,i,startid);
  count = endid-startid+1;
  for(j=startid;j<=endid;j++)
   {
   seg2list[j] = i+(j-startid)*1.0/(count);
   }
  startid = endid+1;
  }
 }


for(i=0;i<chrnum;i++)
 {
 shiftvalue = seg1list[chrbid1[i]]-seg2list[chrbid2[i]];
// printf("%d %d %d %lf %lf %lf\n",i,chrbid1[i],chrbid2[i],seg1list[chrbid1[i]],seg2list[chrbid2[i]],shiftvalue);
 for(j=0;j<totallen2;j++)
  {
  if(chrassign2[assign2[j]] == i)
   seg2list[j] += shiftvalue;
  }
 }



//for(i=0;i<totallen1;i++)
// printf("seg1 %d %d %lf\n",i,chrassign1[assign1[i]],seg1list[i]);
//
//for(i=0;i<totallen2;i++)
// printf("seg2 %d %d %lf\n",i,chrassign2[assign2[i]],seg2list[i]);


tempid1 = 0;
tempid2 = 0;

for(i=0;i<totallen2;i++)
 {
// printf("check %d %d %d %d\n",i,totallen2,tempid1,totallen1);
 tempid1 = getseglistid(seg1list,totallen1,seg2list[i]);
// printf("test %d %lf %d %lf %d %d\n",i,seg2list[i],tempid1,seg1list[tempid1],chrassign1[assign1[tempid1]],chrassign2[assign2[i]]);
 for(count=0;chrassign1[assign1[tempid1]] != chrassign2[assign2[i]];count++)
  {
  if(count > 10)
   {
   printf("ERROR in the data\n");
   exit(0);
   }
  if(chrassign1[assign1[tempid1]] < chrassign2[assign2[i]])
   tempid1 ++;
  else
   tempid1 --;

//  printf("%d %d\n",count,tempid1);

  }
// if(chrassign1[assign1[tempid1]] != chrassign2[assign2[i]])

// printf("check %d %lf %d %lf %d %d\n",i,seg2list[i],tempid1,seg1list[tempid1],chrassign1[assign1[tempid1]],chrassign2[assign2[i]]);
 if(tempid1 == totallen1-1)
  tempid1 --;
 if(chrassign1[assign1[tempid1]] == chrassign1[assign1[tempid1]+1])
  {
  tempdist = seg2list[i]-seg1list[tempid1];
  range = seg1list[tempid1+1]-seg1list[tempid1];
  for(k=0;k<DIM;k++)
   coor2[i][k] = coor[tempid1][k]*(1-tempdist/range)+coor[tempid1+1][k]*tempdist/range;

// printf("%lf %lf %lf\n",coor[tempid1][0],coor[tempid1][1],coor[tempid1][2]);
// printf("1 %d %d %lf %lf %lf\n",statenum,i,coor2[i][0],coor2[i][1],coor2[i][2]);
  }
 else
  {
//  printf("check %d\n",i);
  tempid1 --;
  tempdist = seg2list[i]-seg1list[tempid1];
  range = seg1list[tempid1+1]-seg1list[tempid1];
  for(k=0;k<DIM;k++)
   coor2[i][k] = coor[tempid1][k]*(1-tempdist/range)+coor[tempid1+1][k]*tempdist/range;
// printf("2 %d %lf %lf %lf\n",i,coor2[i][0],coor2[i][1],coor2[i][2]);
  }

 }






outputfile = openfile(outputfilename,'w');
for(i=0;i<totallen2;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",coor2[i][j]*expratio);
 fprintf(outputfile,"\n");
 }


fclose(outputfile);

for(i=0;i<statenum;i++)
 {
 free(coor[i]);
 }
free(coor);

for(i=0;i<totallen2;i++)
 {
 free(coor2[i]);
 }
free(coor2);
free(assign1);
free(assign2);
free(chrlen1);
free(chrlen2);
free(chrassign1);
free(chrassign2);
free(seg1list);
free(seg2list);

free(chrbid1);
free(chrbid2);





}







