#include <stdio.h>
#include <string.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3
#include <math.h>





void main()
{
char inputfilename[N],input2filename[N],outputfilename[N],assignfilename[N];
int i,j,datanum,*assign,assignlen,count,winsize,pointnum,startid,endid,statenum,bid,eid;
double *data,ave,ave2,std,*data2,std2,corr;
FILE *inputfile,*outputfile,*assignfile,*input2file;

printf("Input filename for the input:\n");
scanf("%s",inputfilename);

printf("Input filename for the second data:\n");
scanf("%s",input2filename);

printf("Input the filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input the value for winsize:\n");
scanf("%d",&winsize);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

statenum = getlinenum(assignfilename);

data = doublearray(statenum);
data2 = doublearray(statenum);
assign = intarray(statenum);

assignfile = openfile(assignfilename,'r');
inputfile = openfile(inputfilename,'r');
input2file = openfile(input2filename,'r');
for(i=0;i<statenum;i++)
 {
 fscanf(assignfile,"%d",&assign[i]);
 fscanf(inputfile,"%lf",&data[i]);
 fscanf(input2file,"%lf",&data2[i]);
 }
fclose(assignfile);
fclose(inputfile);
fclose(input2file);

outputfile = openfile(outputfilename,'w');

for(i=0;i<statenum;i++)
 {
 std = 0;
 std2 = 0;
 corr = 0;
 ave = 0;
 ave2 = 0;

 bid = i;
 count = 0;
 eid = i+winsize;
 if(eid > statenum)
  eid = statenum;

 for(j=bid;j<eid;j++)
  {
  if(assign[j] != -1)
   {
   ave += data[j];
   ave2 += data2[j];
   count ++;
   }
  }

 if(count > 0)
  {
  ave /= count;
  ave2 /= count;
  for(j=bid;j<eid;j++)
   {
   if(assign[j] != -1)
    {
    std += (data[j]-ave)*(data[j]-ave);
    std2 += (data2[j]-ave2)*(data2[j]-ave2);
    corr += (data[j]-ave)*(data2[j]-ave2);
    }
   }

  std = sqrt(std/count);
  std2 = sqrt(std2/count);

  if(std*std2 > 0.00001)
   {
   corr /= std*std2*count;
   }
  }

 fprintf(outputfile,"%d %lf %d\n",i+1,corr,count);
 }

fclose(outputfile);


free(data);
free(assign);
free(data2);


}




