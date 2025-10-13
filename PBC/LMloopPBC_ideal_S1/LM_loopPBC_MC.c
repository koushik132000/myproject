//Program provides lattice moves of 36 chains in 7x7x20 lattice
//Program introduces maxm mc steps and moves are accepted or rejected based on the metropolis algorithm
//Creation Date:27Nov2024
//Modified Date:11Oct2025
	//01-12-2024: Introduced subroutine to determine state of a lattice site
 	//05-12-2024: Open-MPI parallelisation is introduced
	//06-12-2024: Pragma for directions introduced. fpos and npos have 6 direction values
	//06-12-2024: Energy change calculations subroutine introduced
	//07-12-2024: Direction evaluation subroutine introduced
	//07-12-2024: Remove array formulations for terminal bead motion
	//07-12-2024: Introduce overlap energy calculations
	//08-12-2024: Introduced check for number of beads occupying a site
	//09-12-2024: Corrected neighborlist calculation in deltaE subroutine
	//10-12-2024: Incorporating periodic boundary condition (PBC)
	//10-01-2025: Writing file after move increment is done
	//21-02-2025: Updated comments for the program
	//06-02-2025: Updated for variable number of simulation steps and dumping data every mcsteps
 	//11-10-2025: Correction for energy calculations

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ran2.h"

#define nc 36//number of chains
#define bc 20//number of beads per chain
#define lx 7//box size in x-direction
#define ly 7//box size in y-direction
#define lz 20//box size in z-direction
#define lxly 49//number of sites in xy layer
#define ns 980//number of sites 

//Subroutines
//sstate: Determine lattice site status and maps it to chain details
//fmoves: Moves the first bead of the chain
//nmoves: Moves the last bead of the chain
//kmoves: Moves kth internal bead of the chain
//deltaE: Energy change associated with a move of any bead
//fmeval: Evaluates energy associated with move of first bead
//lmeval: Evaluates energy associated with move of last bead
//metrop: Metropolis algorithm
//accmov: Accepts moves and reconfigures the chain

struct pos{
	int x[bc];
	int y[bc];
	int z[bc]; 
};

struct vec{
	int ex;
	int ey;
	int ez;
};

struct site{
	int sx[ns];
	int sy[ns];
	int sz[ns];
};

//Chain beads position for nc chains
struct pos beads[nc];

//Lattice configuration and state for 7x7x20 sites
struct site lcs;

//Defining unit vectors for the directions in the cubic lattice
struct vec dir[6]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
struct vec fpos;
struct vec npos;
struct vec kpos;

int nb=bc-1;//index of the last bead
int pb=bc-2;//index of the penultimate bead
int cb=bc-3;//index of bead connected to penultimate bead

int scalc;//Calculated site index
int pcalc;//Previous site index
int fcalc,lcalc;//First and last bead site index
double dEf,dEl,dE;//Energy difference

double Ea=0.0;//Bead-Wall interaction energy
double Eb=0.00;//Bead-Bead interaction energy
double Ex=200.0;//Bead overlap energy (Exclude Volume)

long seed;//seed for the ran2 generator
	
//Determines lattice site details
void sstate()
{
	int i,j,k;
	int xval,yval,zval;

	//Initialize lattice site status and details
	for (k=0;k<ns;k++){
		lcs.sx[k]=0;
		lcs.sy[k]=0;
		lcs.sz[k]=0;
	}

	//Update occupation status and lattice site details
	for (i=0;i<nc;i++){
		for (j=0;j<bc;j++){
			xval=beads[i].x[j];
			yval=beads[i].y[j];
			zval=beads[i].z[j];			
			k=xval+lx*yval+lxly*zval;
			lcs.sx[k]=i;
			lcs.sy[k]=j;
			lcs.sz[k]+=1;
		}
	}
	
}

//First bead of a chain is moved
void fmoves(int cnum,int r)
{
	int move;
	struct vec fvec;

	//Determing the current orientation of the first segment
	fvec.ex=beads[cnum].x[1]-beads[cnum].x[0];
	fvec.ey=beads[cnum].y[1]-beads[cnum].y[0];
	fvec.ez=beads[cnum].z[1]-beads[cnum].z[0];

	//Determining if the bead should be moved
	if (fvec.ex==dir[r].ex && fvec.ey==dir[r].ey && fvec.ez==dir[r].ez)  
		move=0;
	else
		move=1;
		
	//Moving the first bead
	if (move==1){
		fpos.ex=beads[cnum].x[1]-dir[r].ex;
		fpos.ey=beads[cnum].y[1]-dir[r].ey;
		fpos.ez=beads[cnum].z[1]-dir[r].ez;	
	}
	else{
		fpos.ex=beads[cnum].x[0];
		fpos.ey=beads[cnum].y[0];
		fpos.ez=beads[cnum].z[0];
	}
	
	//Periodic boundary condition
	if(fpos.ex<0)
		fpos.ex=lx-1;
	if(fpos.ey<0)
		fpos.ey=ly-1;
	if(fpos.ez<0)
		fpos.ez=lz-1;
	if(fpos.ex>lx-1)
		fpos.ex=0;
	if(fpos.ey>ly-1)
		fpos.ey=0;
	if(fpos.ez>lz-1)
		fpos.ez=0;

}

//Last bead of a chain is moved
void nmoves(int cnum,int r)
{
	int move;
	struct vec nvec;

	//Determing the current orientation of the last segment
	nvec.ex=beads[cnum].x[nb]-beads[cnum].x[pb];
	nvec.ey=beads[cnum].y[nb]-beads[cnum].y[pb];
	nvec.ez=beads[cnum].z[nb]-beads[cnum].z[pb];

	//Determining if the bead should be moved
	if (nvec.ex==dir[r].ex && nvec.ey==dir[r].ey && nvec.ez==dir[r].ez)  
		move=0;
	else
		move=1;
	
	//Moving the last bead
	if (move==1){
		npos.ex=beads[cnum].x[pb]+dir[r].ex;
		npos.ey=beads[cnum].y[pb]+dir[r].ey;
		npos.ez=beads[cnum].z[pb]+dir[r].ez;	
	}
	else{	
		npos.ex=beads[cnum].x[nb];
		npos.ey=beads[cnum].y[nb];
		npos.ez=beads[cnum].z[nb];	
	}
	
	//Periodic boundary condition
	if(npos.ex<0)
		npos.ex=lx-1;
	if(npos.ey<0)
		npos.ey=ly-1;
	if(npos.ez<0)
		npos.ez=lz-1;
	if(npos.ex>lx-1)
		npos.ex=0;
	if(npos.ey>ly-1)
		npos.ey=0;
	if(npos.ez>lz-1)
		npos.ez=0;

}

//kth non-terminal bead of a chain is moved
void kmoves(int cnum, int k)
{
	int dota,move;
	struct vec vec1;
	struct vec vec2;

	//Determing the current orientation of the kth segement
	vec1.ex=beads[cnum].x[k+1]-beads[cnum].x[k];
	vec1.ey=beads[cnum].y[k+1]-beads[cnum].y[k];
	vec1.ez=beads[cnum].z[k+1]-beads[cnum].z[k];
	
	//Determing the current orientation of the (k-1)th segement
	vec2.ex=beads[cnum].x[k]-beads[cnum].x[k-1];
	vec2.ey=beads[cnum].y[k]-beads[cnum].y[k-1];
	vec2.ez=beads[cnum].z[k]-beads[cnum].z[k-1];

	//Determining if the bead should be moved
	dota=vec1.ex*vec2.ex+vec1.ey*vec2.ey+vec1.ez*vec2.ez;
	if (dota!=0)  
		move=0;
	else
		move=1;
	
	//Moving the chosen bead
	if (move==1){
		kpos.ex=beads[cnum].x[k-1]+vec1.ex;
		kpos.ey=beads[cnum].y[k-1]+vec1.ey;
		kpos.ez=beads[cnum].z[k-1]+vec1.ez;	
	}
	else{
		kpos.ex=beads[cnum].x[k];
		kpos.ey=beads[cnum].y[k];
		kpos.ez=beads[cnum].z[k];	
	}
	
	//Periodic boundary condition
	if(kpos.ex<0)
		kpos.ex=lx-1;
	if(kpos.ey<0)
		kpos.ey=ly-1;
	if(kpos.ez<0)
		kpos.ez=lz-1;
	if(kpos.ex>lx-1)
		kpos.ex=0;
	if(kpos.ey>ly-1)
		kpos.ey=0;
	if(kpos.ez>lz-1)
		kpos.ez=0;

}

//Energy change calculation
double deltaE(int olds,int news)
{
	int olist[6],nlist[6];
	int oldn,newn;
	int xy,xz,yz;
	int i;
	double Ediff;
	double overlap_penalty = 0.0;

	/*if(lcs.sz[news]>=1){
		Ediff = lcs.sz[news]*Ex;
		//printf("overlap\n");
		return Ediff;
	}*/
	
	//olds and news are the site indices of a lattice site
	//xy plane neighbors, xz plane neighbors and yz plane neighbors need to be identified
	//xy bottom & top, xz front & back, yz left & right need special considerations
	//xy has 0-19 layers, xz has 0-6 layers, yz has 0-6 layers
	xy=floor(olds/lxly);
	xz=floor((olds-xy*lxly)/ly);
	yz=olds-xy*lxly-xz*ly;
	
	//Neighbors of the old site
	//if((olds+1)%lx!=0)
	//Detemining x neighbors
	if(yz!=lx-1)//not right boundary
		olist[0]=olds+1;
	else
		olist[0]=olds-(lx-1);
	//if(olds%lx!=0)
	if(yz!=0)//not left boundary
		olist[1]=olds-1;
	else
		olist[1]=olds+(lx-1);
	//if((olds+1)%lxly!=0)
	//Detemining y neighbors
	if(xz!=ly-1)//not back boundary
		olist[2]=olds+lx;
	else
		olist[2]=yz*lx;
	//if(olds%lxly!=0)
	if(xz!=0)//not front boundary
		olist[3]=olds-lx;
	else
		olist[3]=olds+lx*(ly-1);
	//if(olds<931)
	//Determining z neighbors
	if(xy!=lz-1)//not top boundary
		olist[4]=olds+lxly;
	else
		olist[4]=olds-lxly*(lz-1);
	//if(olds>(lxly-1))
	if(xy!=0)//not bottom boundary
		olist[5]=olds-lxly;
	else
		olist[5]=olds+lxly*(lz-1);
	
	//xy plane neighbors, xz plane neighbors and yz plane neighbors need to be identified
	//xy bottom & top, xz front & back, yz left & right need special considerations
	//xy has 0-19 layers, xz has 0-6 layers, yz has 0-6 layers
	xy=floor(news/lxly);
	xz=floor((news-xy*lxly)/ly);
	yz=news-xy*lxly-xz*ly;
	
	//Neighbors of the new site
	//if((news+1)%7!=0)
	//Determining x neighbors
	if(yz!=lx-1)//not right boundary
		nlist[0]=news+1;
	else
		nlist[0]=news-(lx-1);
	//if(news%7!=0)
	if(yz!=0)//not left boundary
		nlist[1]=news-1;
	else
		nlist[1]=news+(lx-1);
	//if((news+1)%49!=0)
	//Determining y neighbors
	if(xz!=lx-1)//not back boundary
		nlist[2]=news+lx;
	else
		nlist[2]=yz*lx;
	//if(news%49!=0)
	if(xz!=0)//not front boundary
		nlist[3]=news-lx;
	else
		nlist[3]=news+lx*(ly-1);
	//if(news<931)
	//Determining z neighbors
	if(xy!=lz-1)//not top boundary	
		nlist[4]=news+lxly;
	else
		nlist[4]=news-lxly*(lz-1);
	//if(news>48)
	if(xy!=0)//not bottom boundary	
		nlist[5]=news-lxly;
	else
		nlist[5]=news+lxly*(lz-1);

	//Number of neighbors in the old and new state
	oldn=0;
	newn=0;
	for (i=0;i<6;i++){
		if(olist[i]>=0 && olist[i]<980){
			oldn+=lcs.sz[olist[i]];
		}
		if(nlist[i]>=0 && nlist[i]<980){
			newn+=lcs.sz[nlist[i]];
		}	
	}
	oldn-=1;
	newn+=1;
	if (Ex != 0.0 && lcs.sz[news] >= 1) {
		overlap_penalty = lcs.sz[news] * Ex;
	}
	//Energy difference
	Ediff=(newn-oldn)*Eb + overlap_penalty;
	
	return Ediff;
}

//Evaluates the terminal moves for first bead in a specific direction	
void fmeval(int cnum,int dindex,int fcalc)
{
	//Energy calculations for moving the chosen terminal bead in dindex
	fmoves(cnum,dindex);
	scalc=fpos.ex+lx*fpos.ey+lxly*fpos.ez;
	dEf=deltaE(fcalc,scalc);
}	

//Evaluates the terminal move for the last bead in aspecific direction
void lmeval(int cnum,int dindex,int lcalc)
{
	nmoves(cnum,dindex);
	scalc=npos.ex+lx*npos.ey+lxly*npos.ez;
	dEl=deltaE(lcalc,scalc);
}

//Metropolis algorithm implementation
int metrop(double delE)
{
	int acc;
	double rij;
	double pij;

	if (delE<=0){
		pij=1.0;
		acc=1;
	}
	else{
		pij=exp(-delE);
		rij=ran2(&seed);
		if (pij<rij)
			acc=0;
		else
			acc=1;
	}
	return acc;
}

//Accepts move and reconfigures the chain
void accmov(int cnum,int bnum,int pcalc,int scalc,struct vec mpos)
{
	lcs.sx[pcalc]=0;
	lcs.sy[pcalc]=0;
	lcs.sz[pcalc]-=1;
	beads[cnum].x[bnum]=mpos.ex;
	beads[cnum].y[bnum]=mpos.ey;
	beads[cnum].z[bnum]=mpos.ez;
	lcs.sx[scalc]=cnum;
	lcs.sy[scalc]=bnum;
	lcs.sz[scalc]+=1;
}

//Main simulation
int main()
{
	FILE *fptr;
	FILE *fptw;

	char x,y,z;
	int ccount=0;
	int bcount=0;
	int read=0;
	int sindex;//Site index
	int findex,lindex;//chosen move directions
	int cnum;//chosen chain
	int bnum;//chosen bead
	int i,j;
	int fn;//file number
	int sucm,totm;//number of mc moves
	int maxm=200E6;
	int mcstep=5E6;

	char chain_file[15];
	char file_number[3];

	//Reading data from the input file
	fptr=fopen("chains.csv","r");
	if (fptr== NULL)
	{
		printf("Error reading file\n");
		return 1;
	}
	do{
		if (bcount==0 && ccount==0)
		{
			read=fscanf(fptr,"%c,%c,%c\n",&x,&y,&z);//File header read
		}	
		//Data reading
		read=fscanf(fptr,"%d,%d,%d\n",&beads[ccount].x[bcount],&beads[ccount].y[bcount],&beads[ccount].z[bcount]);
		if (read==3)
			bcount = bcount+1;
		else
			printf("read error %d\n",read);
		if (bcount==bc)
		{
			ccount=ccount+1;
			bcount=0;
		}
			
	}while (!feof(fptr));
	fclose(fptr);

	/*srand(time(NULL));
	seed=(long)time(NULL);*/
	seed = -12345;
	if (seed>=0)
		seed=-1-seed;
	for (i=0;i<100;i++)
		printf("%f %f\n",(double)rand()/RAND_MAX,ran2(&seed));

	//Determining lattice state in the simulation box
	sstate();
	for (sindex=0;sindex<ns;sindex++){
		if (lcs.sz[sindex]==1)
			printf("site: %d %d %d %d\n",sindex,lcs.sx[sindex],lcs.sy[sindex],lcs.sz[sindex]);
		else
			printf("site: %d not occupied\n",sindex);
	}
	fn=0;
	sucm=0;
	totm=0;
	//Write initial configuration to a file 
        sprintf(chain_file,"mchains%d.csv",fn);
        fptw=fopen(chain_file,"wb");
        fprintf(fptw,"x,y,z\n");
        for (i=0; i<nc; i++){
                for (j=0; j<bc; j++){
                        fprintf(fptw,"%d,%d,%d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
                }
        }
        fclose(fptw);
	do{			
		//Randomly choose a chain	
		cnum =(int)(ran2(&seed)*nc);
		//Move the first bead
		dEf=0.0;
		fcalc=beads[cnum].x[0]+lx*beads[cnum].y[0]+lxly*beads[cnum].z[0];
		findex=(int)(ran2(&seed)*6);
		fmeval(cnum,findex,fcalc);
		if(metrop(dEf)==1){	
			scalc=fpos.ex+lx*fpos.ey+lxly*fpos.ez;
			accmov(cnum,0,fcalc,scalc,fpos);
			sucm+=1;
		}
		//Move last bead
		dEl=0.0;
		lcalc=beads[cnum].x[nb]+lx*beads[cnum].y[nb]+lxly*beads[cnum].z[nb];
		lindex=(int)(ran2(&seed)*6);
		lmeval(cnum,lindex,lcalc);
		if (metrop(dEl)==1){
			scalc=npos.ex+lx*npos.ey+lxly*npos.ez;	
			accmov(cnum,nb,lcalc,scalc,npos);
			sucm+=1;
		}
		//Move penultimate bead
		dE=0.0;	
		pcalc=beads[cnum].x[pb]+lx*beads[cnum].y[pb]+lxly*beads[cnum].z[pb];	
		kmoves(cnum,pb);
		scalc=kpos.ex+lx*kpos.ey+lxly*kpos.ez;
		dE=deltaE(pcalc,scalc);	
		if (metrop(dE)==1){
			accmov(cnum,pb,pcalc,scalc,kpos);
			sucm+=1;
		}
		//Move second bead
		dE=0.0;
		pcalc=beads[cnum].x[1]+lx*beads[cnum].y[1]+lxly*beads[cnum].z[1];
		kmoves(cnum,1);
		scalc=kpos.ex+lx*kpos.ey+lxly*kpos.ez;
		dE=deltaE(pcalc,scalc);
		if (metrop(dE)==1){	
			accmov(cnum,1,pcalc,scalc,kpos);
			sucm+=1;
		}
		//Move random internal bead of chain
		do {
			bnum=(int)(ran2(&seed)*bc);
		}while(bnum==0||bnum==19);
		dE=0.0;
		pcalc=beads[cnum].x[bnum]+lx*beads[cnum].y[bnum]+lxly*beads[cnum].z[bnum];
		kmoves(cnum,bnum);
		scalc=kpos.ex+lx*kpos.ey+lxly*kpos.ez;
		dE=deltaE(pcalc,scalc);
		if (metrop(dE)==1){	
			accmov(cnum,bnum,pcalc,scalc,kpos);
			sucm+=1;
		}
		totm+=1;
		//Write to file
		if(totm%mcstep==0){
			printf("I am here\n");
			fn+=1;
			sprintf(chain_file,"mchains%d.csv",fn);
			fptw=fopen(chain_file,"wb");
			fprintf(fptw,"x,y,z\n");
			for (i=0; i<nc; i++){
				for (j=0; j<bc; j++){
					fprintf(fptw,"%d,%d,%d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
				}
			}
			fclose(fptw);
		}
	}while(totm<maxm);

	printf("Total number of moves %d and successful of moves %d\n",totm,sucm);
	//Check if a site has more than one bead
	sstate();
	for (sindex=0;sindex<ns;sindex++){
		if (lcs.sz[sindex]>1)
			printf("site: %d %d %d %d\n",sindex,lcs.sx[sindex],lcs.sy[sindex],lcs.sz[sindex]);
	}

	fptw = fopen("mchainsf.csv","w+");
	fprintf(fptw,"x,y,z\n");
	for (i=0; i<nc; i++){
		for (j=0; j<bc; j++){
			fprintf(fptw,"%d,%d,%d\n",beads[i].x[j],beads[i].y[j],beads[i].z[j]);
		}
	}
	fclose(fptw);

	return 0;
}
