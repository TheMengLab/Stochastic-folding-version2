#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3

void main()
{
char microfilename[N],mapfilename[N],macrofilename[N];
int i,j,**map,micronum,linenum;
FILE *microfile,*mapfile,*macrofile;

printf("Input filename for original micro assignment:\n");
scanf("%s",microfilename);

printf("Input filename for micro2macro map(only one column):\n");
scanf("%s",mapfilename);

printf("Input filename for output:\n");
scanf("%s",macrofilename);

micronum = getlinenum(mapfilename);
linenum = getlinenum(microfilename);

map = intmatrix(micronum,2);

mapfile = openfile(mapfilename,'r');
for(i=0;i<micronum;i++)
 for(j=0;j<2;j++)
  fscanf(mapfile,"%d",&map[i][j]);
fclose(mapfile);


microfile = openfile(microfilename,'r');
macrofile = openfile(macrofilename,'w');

for(i=0;i<linenum;i++)
 {
 fscanf(microfile,"%d",&j);
 printf("%d\n",j);
 fprintf(macrofile,"%d %d\n",map[j][0],map[j][1]);
// printf("after %d\n",j);
 }
fclose(microfile);
fclose(macrofile);

for(i=0;i<micronum;i++)
 free(map[i]);

free(map);



}




