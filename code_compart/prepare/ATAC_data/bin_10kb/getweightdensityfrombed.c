#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>





void adddensity(double *density,double binsize,int statenum,int startpoint,int endpoint)
 {
 int i,j,startbin,endbin;

 startbin = (int)(startpoint/binsize);
 endbin = (int)(endpoint/binsize);

 if(startbin == endbin)
  {
  density[startbin] += (endpoint-startpoint+1)/binsize;
  }
 else
  {
  density[startbin] += ((startbin+1)*binsize-startpoint)/binsize;
  for(i=startbin+1;i<=endbin-1;i++)
   density[i] += 1;
  density[endbin] += (endpoint-endbin*binsize)/binsize;
  }

 }






void addweightdensity(double *density,double binsize,int statenum,int startpoint,int endpoint,double weight)
 {
 int i,j,startbin,endbin;

 startbin = (int)(startpoint/binsize);
 endbin = (int)(endpoint/binsize);

 if(startbin == endbin)
  {
  density[startbin] += weight;
  }
 else
  {
  density[startbin] += ((startbin+1)*binsize-startpoint)*weight/(double)(endpoint-startpoint);
  for(i=startbin+1;i<endbin;i++)
   {
   density[i] += binsize*weight/(double)(endpoint-startpoint);
   }
  density[endbin] += (endpoint-endbin*binsize)*weight/(endpoint-startpoint);
//  for(i=startbin+1;i<=endbin-1;i++)
//   density[i] += 1;
//  density[endbin] += (endpoint-endbin*binsize)/binsize;
  }

 }




void addweightdensity_bin(double *density,double binsize,int statenum,int currentid,int *startlist,int *endlist,int linenum,double *weight,double shiftwinsize)
 {
 int i,j;
 double startid,endid;

 startid = currentid*binsize-shiftwinsize;
 endid = (currentid+1)*binsize+shiftwinsize;

 for(i=0;i<linenum;i++)
  {
  if((startid > endlist[i])||(endid < startlist[i]))
   continue;
  else
   {
   if(((startid>=startlist[i])&&(startid<endlist[i]))&&(endid>=endlist[i]))
    {
 //   density[currentid] += weight[i]*(endlist[i]-startid)/(endlist[i]-startlist[i]);
    density[currentid] += weight[i]*(endlist[i]-startid)/(endid-startid);
 //   printf("%d %d %d %d %lf %lf\n",currentid,startlist[i],endlist[i],endlist[i]-startlist[i],weight[i],weight[i]*(endlist[i]-startid)/(endid-startid));
    }
   else if(((startid>=startlist[i])&&(startid<endlist[i]))&&(endid<endlist[i]))
    {
    density[currentid] += weight[i];
 //   printf("%d %d %d %d %lf %lf\n",currentid,startlist[i],endlist[i],endlist[i]-startlist[i],weight[i],weight[i]*(endlist[i]-startid)/(endid-startid));
    }
   else if(((endid>=startlist[i])&&(endid<endlist[i]))&&(startid<startlist[i]))
    {
 //   density[currentid] += weight[i]*(endid-startlist[i])/(endlist[i]-startlist[i]);
    density[currentid] += weight[i]*(endid-startlist[i])/(endid-startid);
 //   printf("%d %d %d %d %lf %lf\n",currentid,startlist[i],endlist[i],endlist[i]-startlist[i],weight[i],weight[i]*(endlist[i]-startid)/(endid-startid));
    }
   else if((startid < startlist[i])&&(endid >= endlist[i]))
    {
    density[currentid] += weight[i]*(endlist[i]-startlist[i])/(endid-startid);
    }
   }
  }
 }



void main()
{
char inputfilename[N],outputfilename[N];
int i,j,k,linenum,start,end,statenum,*startlist,*endlist,max;
double binsize,*density,*weight,shiftwinsize;
FILE *inputfile,*outputfile;

printf("Input the filename for the range:\n");
scanf("%s",inputfilename);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

binsize = 10000.0;
shiftwinsize = 0.0;
//shiftwinsize = 175000.0;
linenum = getlinenum(inputfilename);

startlist = intarray(linenum);
endlist = intarray(linenum);
weight = doublearray(linenum);

max = 0;
inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&startlist[i]);
 fscanf(inputfile,"%d",&endlist[i]);
 fscanf(inputfile,"%lf",&weight[i]);
 if(max < endlist[i])
  max = endlist[i];
 }
fclose(inputfile);

statenum = (int)(max/binsize) +1;

density = doublearray(statenum);

for(i=0;i<statenum;i++)
 density[i] = 0;

for(i=0;i<statenum;i++)
 {
 addweightdensity_bin(density,binsize,statenum,i,startlist,endlist,linenum,weight,shiftwinsize);
 printf("%d\n",i);
 }

//for(i=0;i<linenum;i++)
// {
//// adddensity(density,binsize,statenum,startlist[i],endlist[i]);
// addweightdensity(density,binsize,statenum,startlist[i],endlist[i],weight[i]);
// }

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 fprintf(outputfile,"%lf\n",density[i]);
fclose(outputfile);

free(startlist);
free(endlist);
free(density);
free(weight);






}







