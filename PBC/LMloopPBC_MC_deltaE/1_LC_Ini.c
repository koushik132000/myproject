//Program provides initial configuration of 36 chains in 7x7x20 lattice
//Creation Date:30Oct2024
//Last Modified:10Nov2024
//Modifications:
	//07-11-2024:Introduced 3 subroutines for arranging chains on the xy-plane 
	//(1) along diagonals
	//(2) aligned along x direction
	//(3) aligned along y direction
	//10-11-2024:Introduced subroutine for arranging chains randomly

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#define nc 36//number of chains
#define bc 20//number of beads per chain
#define lx 7//box size in x-direction
#define ly 7//box size in y-direction
#define lz 20//box size in z-direction

//Subroutines
//diagonal: Arranges the chains along the diagonal of xy-plane
//xarrange: Arrange chains along x-direction
//yarrange: Arrange chains along y-direction
//randomize:Arrange chains randomly in the xy-plane
 
struct pos{
	int x[bc];
	int y[bc];
	int z[bc]; 
};

struct pos beads[nc];

//Subrotuine for arranging rigid chains along the diaganol of xy-plane
void diagonal()
{
	int i,j;
	int lsx,lsy;

	i=0;//chain number
	lsx=0;
	lsy=0;
	
	while (i<nc)
	{	

		printf("Chains number %d bead positions \n",i+1);
		for(j=0;j<bc;j++)
		{
			beads[i].x[j]=lsx;
			beads[i].y[j]=lsy;
			beads[i].z[j]=j;
			printf("%d %d %d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
		i=i+1;
		//Diagonal
		if (i<lx)
		{
			lsx=i%lx;
			lsy=i%ly;
		}
		//Diagonal-1
		if (i>=lx && i<2*lx-1)
		{
			lsx=i%lx+1;
			lsy=i%ly;
		}
		//Diagonal+1
		if (i>=2*lx-1 && i<3*lx-2)
		{
			lsx=(i+1)%lx;
			lsy=(i+1)%ly+1;
		}
		//Diagonal-2
		if (i>=3*lx-2 && i<4*lx-4)
		{
			lsx=(i+2)%lx+2;
			lsy=(i+2)%ly;
		}
		//Diagonal+2
		if (i>=4*lx-4 && i<5*lx-6)
		{
			lsx=(i+4)%lx;
			lsy=(i+4)%ly+2;
		}
		//Diagonal-3
		if (i>=5*lx-6 && i<6*lx-9)
		{
			lsx=(i+6)%lx+3;
			lsy=(i+6)%ly;
		}
		//Diagonal+3
		if (i>=6*lx-9 && i<7*lx-12)
		{
			lsx=(i+9)%lx;
			lsy=(i+9)%ly+3;
		}
	}
}

//Subrotuine for arranging rigid chains along the x-axis of xy-plane
void xarrange()
{
	int i,j;
	int lsx,lsy;

	i=0;//chain number
	lsx=0;
	lsy=0;
	
	while (i<nc)
	{	
		printf("Chains number %d bead positions \n",i+1);
		for(j=0;j<bc;j++)
		{
			beads[i].x[j]=lsx;
			beads[i].y[j]=lsy;
			beads[i].z[j]=j;
			printf("%d %d %d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
		i=i+1;
		if(i%lx==0)
		{
			lsx=lsx+1;
		}
		lsy=i%lx;
	}
}

//Subrotuine for arranging rigid chains along the y-axis of xy-plane
void yarrange()
{
	int i,j;
	int lsx,lsy;

	i=0;//chain number
	lsx=0;
	lsy=0;
	
	while (i<nc)
	{	
		printf("Chains number %d bead positions \n",i+1);
		for(j=0;j<bc;j++)
		{
			beads[i].x[j]=lsx;
			beads[i].y[j]=lsy;
			beads[i].z[j]=j;
			printf("%d %d %d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
		i=i+1;
		lsx=i%lx;
		if(i%ly==0)
		{
			lsy=lsy+1;
		}
	}
}

//Subrotuine for arranging rigid chains randomly on the xy-plane
void randomize()
{
	int i,j,k;
	int valuex[nc],valuey[nc];

	srand(0);//seeding the random number
	i=0;//chain number initialization
	while(i < nc)
	{
		valuex[i]=rand()%lx;
		valuey[i]=rand()%ly;
		printf("Sites (x,y)=(%d,%d)\n",valuex[i],valuey[i]);
		for (k=0;k<i;k++)
		{ 
			if (valuex[i]== valuex[k] && valuey[i]==valuey[k])
			{
				printf("Repeated sites is true reset counter\n");
				i=i-1;
			}
		}
		i=i+1;
	}

	i=0;	
	while (i<nc)
	{	
		printf("Chains number %d bead positions \n",i+1);
		for(j=0;j<bc;j++)
		{
			beads[i].x[j]=valuex[i];
			beads[i].y[j]=valuey[i];
			beads[i].z[j]=j;
			printf("%d %d %d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
		i=i+1;
	}
}


//Writes chain configuration to csv file
int main()
{
	int i,j;
	FILE *fpt;

	
	fpt = fopen("chains.csv","w+");
	fprintf(fpt,"x,y,z\n");
	
	//diagonal();
	xarrange();
	//yarrange();
	//randomize();

	for (i=0; i<nc; i++)
	{
		for (j=0; j<bc; j++)
		{	
			fprintf(fpt,"%d,%d,%d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
	}

	fclose(fpt);

	return 0;
}
