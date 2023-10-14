#include <stdio.h>
#include <string.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>





void main()
{
char inputfilename[N],outputfilename[N],assignfilename[N];
int i,j,datanum,*assign,assignlen,count;
double *data,average,std;
FILE *inputfile,*outputfile,*assignfile;

printf("Input filename for the input:\n");
scanf("%s",inputfilename);

printf("Input the filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

datanum = getfloatnum(inputfilename);
assignlen = getintnum(assignfilename);

data = doublearray(datanum);
assign = intarray(assignlen);

assignfile = openfile(assignfilename,'r');
for(i=0;i<assignlen;i++)
 fscanf(assignfile,"%d",&assign[i]);
fclose(assignfile);



inputfile = openfile(inputfilename,'r');
average = 0;
count = 0;
for(i=0;i<datanum;i++)
 {
 fscanf(inputfile,"%lf",&data[i]);
 if(i < assignlen)
  {
  if(assign[i] >= 0)
   {
   average += data[i];
   count ++;
   }
  }
 }

fclose(inputfile);




average /= count;

std = 0;
if(datanum > 1)
 {
 for(i=0;i<datanum;i++)
  if(i<assignlen)
   if(assign[i] >= 0)
    std += (data[i]-average)*(data[i]-average);
 std = sqrt(std/(count-1));
 }

outputfile = openfile(outputfilename,'w');
fprintf(outputfile,"%lf   %lf\n",average,std);

fclose(outputfile);
free(data);
free(assign);


}




