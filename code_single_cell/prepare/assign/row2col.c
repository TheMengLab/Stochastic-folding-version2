#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200

void main()
{
char cptfilename[N],outputfilename[N];
int i,j,k,cptnum,*cptid;
FILE *cptfile,*outputfile;

printf("Input the filename for the cpt file:\n");
scanf("%s",cptfilename);

printf("Input the filename for the output filename:\n");
scanf("%s",outputfilename);

cptnum = getintnum(cptfilename);

cptid = intarray(cptnum);

cptfile = openfile(cptfilename,'r');
outputfile = openfile(outputfilename,'w');

for(i=0;i<cptnum;i++)
 fscanf(cptfile,"%d",&cptid[i]);

for(i=0;i<cptnum;i++)
 fprintf(outputfile,"%d\n",cptid[i]);

//free(hnum);
free(cptid);

//fclose(hfile);
fclose(cptfile);
fclose(outputfile);

}

