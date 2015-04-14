
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "stringlib.h"

// Declare Subroutines
void read_cfg_file(int *, float *, float *, char *, char *, int *, int *, float *);
void write_xyz_step(float **, int, int, float, FILE *);
void init_positions(float **, int, float);
float total_pair_energy(float **, int, float);
float lennard_jones(float *, float *, float);

static float kB = 1.9872041E-3; // actually R in units of kcal/mol/K
static float eps = 0.210849;    // units of kcal/mol
static float sigma = 3.345;  // units of angstroms
static float sigma6 = 1400.80193382; // units of angstroms^6

//Main Program
int main() {

	int nAtoms;     // Number of atoms
	float temp;     // temperature
	float box;      // cubic box size
	int nIter;      // number of MC iterations
	int deltaWrite; // how often to write coordinates and log info in MC
	float deltaX;   // how big to make the translation attempt

	char trajFileName[1024];  // output trajectory file name
	char logFileName[1024];   // log file name

	float **coord;  // coordinates of particles

	int i, j, k;    // genereic indeces
	int iter;       // MC iteration
	int atom;       // MC selected atom
	float energy;   // energy of the system
	float delta[3]; // added position
	float kBT;
	float newEnergy;
	float deltaE;
	int acceptedMoves;
		
	FILE *xyzOut;

	// read config data from standard in
	read_cfg_file(&nAtoms, &temp, &box, trajFileName, logFileName, &nIter, &deltaWrite, &deltaX);
	kBT = kB*temp;
	printf("kB*T=%f\n",kBT);

	// allocate coordinate array
	coord = (float**) malloc(nAtoms*sizeof(float*));
	for (i=0;i<nAtoms;i++) {
		coord[i] = (float*) malloc(3*sizeof(float));
	}
	// initialize particle positions
	init_positions(coord,nAtoms,box);

	// Compute energy of system
	energy = total_pair_energy(coord,nAtoms,box);

	xyzOut = fopen(trajFileName,"w");

	// Perform MC loop
	for(iter=0;iter<nIter;iter++) {
		if (iter%deltaWrite==0) {
			printf("Step: %10d Energy:%50.5f\n", iter, energy);
			// write positions
			write_xyz_step(coord,nAtoms,iter, box, xyzOut);
		}

		// randomly choose a particle to move
		atom = rand()%nAtoms;

		// compute random translation
		for (i=0;i<3;i++) {
			delta[i] = deltaX*(rand()/((float) RAND_MAX)-0.5);
			coord[atom][i] += delta[i];
		}
		// compute new energy
		newEnergy = total_pair_energy(coord,nAtoms,box);

		deltaE = newEnergy-energy;
		if (exp(-deltaE/kBT)> (rand()/((float) RAND_MAX))) {
			energy = newEnergy;
			acceptedMoves++;
			// check to see if we need to wrap
			for(i=0;i<3;i++) {
				if (coord[atom][i] > box) {
					coord[atom][i] -= box;
				} else if (coord[atom][i]<0) {
					coord[atom][i] += box;
				}
			}
		} else {
			for(i=0;i<3;i++) {
				coord[atom][i]-=delta[i];
			}
		}


	}

	fclose(xyzOut);

}

// Subroutines
//

float total_pair_energy(float **coord, int nAtoms, float box) {

	int atom1;
	int atom2;
	float energy;

	energy=0;
	for (atom1=0;atom1<nAtoms-1;atom1++) {

		for (atom2=atom1+1;atom2<nAtoms;atom2++) {

			energy += lennard_jones(coord[atom1],coord[atom2],box);		

		}

	}

	return energy;

}

float lennard_jones(float *pos1, float *pos2, float box) {

	float dist2;
	float temp;
	float energy;
	float dist6;
	int i;

	// compute the distance between the atoms
	dist2 = 0;
	for (i=0;i<3;i++) {
		temp = pos1[i]-pos2[i];
		// check periodic boundaries
		if (temp< -box/2.0) {
			temp += box;
		} else if (temp > box/2.0) {
			temp -= box;
		}
		dist2 += temp*temp;
	}
	dist6 = dist2*dist2*dist2;

	// compute the energy
	energy = 4*eps*(sigma6*sigma6/(dist6*dist6)-sigma6/dist6);

	return energy;	

}

void init_positions(float **coord, int nAtoms, float box) {

	float cbrtf(float x); // cube root function
	int iBoxD;            // integer box dimension
	float fBoxD;          // float box dimension

	int x, y, z;
	float xPos,yPos,zPos;
	int atomCount;

	// determine how many bins to divide the box into
	iBoxD = (int) cbrtf((float) nAtoms);
	if (iBoxD*iBoxD*iBoxD < nAtoms) {
		iBoxD++;
	}
	// determine the size of the bins
	fBoxD = box/( (float) iBoxD);

	// add a particle in each box
	atomCount=0;
	for(x=0;x<iBoxD;x++) {
		xPos = (x+0.5)*fBoxD;
		for (y=0;y<iBoxD;y++) {
			yPos = (y+0.5)*fBoxD;
			for (z=0;z<iBoxD;z++) {
				if (atomCount<nAtoms) {
					zPos = (z+0.5)*fBoxD;
					coord[atomCount][0]=xPos;
					coord[atomCount][1]=yPos;
					coord[atomCount][2]=zPos;
					atomCount++;
				} else {
					break;
				}
			}
			if (atomCount>=nAtoms) {
				break;
			}
		}
		if (atomCount>=nAtoms) {
			break;
		}
	}
}		

void read_cfg_file(int *nAtoms, float *temp, float *box, char *trajFileName, char *logFileName, int *nIter, int *deltaWrite, float *deltaX) {

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;
	float *rCut;

	while (fgets(buffer,1024,stdin) != NULL) {

		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
//		printf ("First word = %s\n",firstWord);
		if (strncmp(firstWord,"nAtoms",6)==0) {
			*nAtoms = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"nIter",5)==0) {
			*nIter = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"deltaWrite",10)==0) {
			*deltaWrite = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"temperature",11)==0) {
			*temp = atof(string_secondword(buffer));
//		} else if (strncmp(firstWord,"rcut",4)==0) {
//			*rCut = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"box",3)==0) {
			*box = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"deltaX",6)==0) {
			*deltaX = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"trajFile",8)==0) {
			strcpy(trajFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"logFile",7)==0) {
			strcpy(logFileName,string_secondword(buffer));
		}
	

	}	
	
	// Print log file info
	printf("Trajectory file: %s\n",trajFileName);
	printf("log file: %s\n",logFileName);
	printf("nAtoms: %d\n",*nAtoms);
	printf("Temperature: %f\n",*temp);
	printf("nIter: %d\n",*nIter);
	printf("deltaWrite: %d\n",*deltaWrite);
	printf("box dimension: %f\n", *box);
	printf("deltaX (MC translation): %f\n", *deltaX);
//	printf("cutoff: %f\n",*rCut);


}


void write_xyz_step(float **coord, int nAtoms, int iter, float box, FILE *xyzOut) {

	int i, j, k;
	int atom;

	fprintf(xyzOut,"%d\n",nAtoms);
	fprintf(xyzOut,"Step %d box %8.3f %8.3f %8.3f\n",iter, box, box, box);
	for (atom=0;atom<nAtoms;atom++) {

		fprintf(xyzOut,"Ar %12.6f%12.6f%12.6f\n",coord[atom][0],coord[atom][1],coord[atom][2]);

	}


}

