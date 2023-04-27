#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "string.h"
#include "math.h"


double randZerotoOne()
{
	double temp = rand();
	double t1 = temp / (double) RAND_MAX;	
	return t1;
}

//Returns approximate value of e^x using sum of first n terms of Taylor Series
float exponential(int n, float x)
{
	float sum = 1.0f; // initialize sum of series
	int i;
	for (i = n - 1; i > 0; --i )
		sum = 1 + x * sum / i;

	return sum;
}

void cellModelLoop(int maxTimeSteps, double sens, int Z, int Zrevert, int ChemoWeeks, int ChemoDays, int ChemoDoses, int ThalfPerDose[7], int ChemoDay[7], int ChemoHour[7], double Cmax, int RadWeeks, int radDoses, int radDoseA[10], int radHour[10], int radDay[10], int *stemCount, int * tbCount, int *cellAge)
{
	int rank, radDose;
	int typeCounter;
	double maxDist;

	srand(1); /* initialize rand(), hardcoded for validation. Change to driven from the time step for full data collection. */	


	/******Variable Declarations*******/
	int DT = 2; //time steps per minute
	char XYZname[35];
	int ratioTbtoS = 20; //ratio of tumor bulk to stem cells before therapy
	double tmp;
	int tau = 3; //number of cell diameters away from vessel wall that stem cells can exit
	double growthRate = .0004;
	int RvesselCells=2; 
	int Rvessel = 4; //radius of the blood vessel
	int R=2; //radius of glial cells
	int k = 5000; //spring constant
	double rTB = .0038; //from Cell paper, proliferation rate of tumor cells after exiting quiescence
	double rS = .0008; //from Cell paper, proliferation rate of stem cells after exiting quiescence

	int timeDepGamma = 0; 
	double mu = 3.25*DT*60; 
	double sigmaS = 1.46*DT*DT*60*60; 
	double lambda = 2*24*60*2; 
	double randomDeathProb=2*24*60*2;
	double aDiff = .0019;
	double ChemoC[7]; 
	double aDeDiff = .45; 
	double randomDiff= 1/aDiff *DT *60;
	double randomDediff=1/aDeDiff *DT *60;
	double Diffusion = .00066;
	double IC50 = .00068; 

	double alpha, alphaS = .00987; 
	double alphaT = .0987; 
	double beta, betaS = 1.14*.00000001; 
	double betaT = 1.14*.0000001; 
	double gamma = .4; 
	double RadRho = .4; 
	double LStem = 36; 
	double LTumorBulk = 24;
	double lambdaStem = .0328; 
	double lambdaTumorBulk =.1; 
	double mTumorBulk = 24*(DT*60); 
	double nuTumorBulk = .54; 

	double hourconv = 1 / (60*DT);
	int chemoDosesGiven = 0; 
	int Thalf = (int) floor(1.1*60*DT*.34); 
	int Tmax =30*DT*.34; 
	double rho;
	int ChemoTime[7];

	ThalfPerDose[0] = -1;
	ThalfPerDose[1] = -1;
	ThalfPerDose[2] = -1;
	ThalfPerDose[3] = -1;
	ThalfPerDose[4] = -1;
	ThalfPerDose[5] = -1;
	ThalfPerDose[6] = -1;

	int giveRad = 0; 
	int RAD_KILL=0;
	int Ch_rm=0;
	int TOT_RAD=0;
	int TOT_DIV = 0;
	int TOT_KILL = 0;

	double ct1= 1; 
	double dt2 = ct1/((double) DT) * ct1/((double) DT);
	int therapyStart = 2800;
	int T = 1; 

	double growth = 1.0057;
	for (T=0; T < maxTimeSteps; T++) {
		tbCount[T] = 0;
		stemCount[T] = 0;
	}

	int simTime = 0; //time elapsed in simulation time, ie DT*T

	double cellRadii[10000];
	double cellsX[10000];
	double cellsY[10000];
	double forcesX[10000];
	double forcesY[10000];
	int cellQ[10000]; 
	int cellDivTime[10000]; 
	double cellConc[10000]; 
	int cellDivRate[10000]; 
	int cellDeathTime[10000]; 
	int cellDiffTime[10000]; 
	int celldeDiffTime[10000]; 
	int cellType[10000]; 
	double cellDist[10000]; 
	int radQ[10000];
	int radQt[10000]; 
	int i, j;
	int numCells = floor( 3.14*2*(Rvessel+2*R)/(2*R)); 

	int max_cycles=1;
	int cycle;
	for (cycle = 0; cycle < max_cycles; cycle++) {
		for (i=0; i < 10000; i++)
		{	
			cellRadii[i] = R;
			cellsX[i] = 0;
			forcesX[i] = 0;
			cellsY[i] = 0;
			forcesY[i] = 0;
			cellQ[i] = 0;
			cellConc[i] = 0;
			cellDivTime[i] = 10000*2;
			cellDeathTime[i] = 10000*2;
			celldeDiffTime[i] = 10000*2;
			cellDiffTime[i] = 10000*2;
			cellDivRate[i] = 0;
			cellType[i] = 0;
			cellDist[i] = 0;
			cellAge[i] = 0;
			radQ[i] = 0;
			radQt[i] = 0;
		}

		int numVesselcells = (int) floor(3.14*2*(Rvessel)/(RvesselCells*2)); //number of vessel cells to start with
		double VesselcellsX[100];
		double VesselcellsY[100];
		for (i=0; i < 100; i++)
		{
			VesselcellsX[i] = 0;
			VesselcellsY[i] = 0;
		}
		for (i=0; i < numVesselcells; i++) {
			VesselcellsX[i] = Rvessel * cos(2*3.14/numVesselcells*i);		
			VesselcellsY[i] = Rvessel * sin(2*3.14/numVesselcells*i);		
		}

		/******Create initial glial cell locations*********/
		numCells = floor( 3.14*2*(Rvessel+2*R)/(2*R)); //number of glial cells to start with
		int countDiv=0, countDeath=0;
		for (i=0; i < numCells; i++) {
			cellsX[i] = (Rvessel+2*R) * cos(2*3.14/numCells*i);
			cellsY[i] = (Rvessel+2*R) * sin(2*3.14/numCells*i);
			cellDivTime[i] = (int) floor(randZerotoOne() * 120); //(lambda));
			cellConc[i] = 0;
			cellDeathTime[i] = (int) floor(randZerotoOne() * (randomDeathProb));
			celldeDiffTime[i] = (int) floor(randZerotoOne() * (randomDediff));
			cellDiffTime[i] =  (int) floor(randZerotoOne() * (randomDiff));
			cellDivRate[i] = (1/rS)*DT*60;
			cellAge[i] = cellDeathTime[i];
			cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])-Rvessel;
			stemCount[0] = stemCount[0] + 1;
			if (cellDeathTime[i] < cellDivTime[i]) countDeath ++;
			else countDiv++;
		}


		/*****Calculate Force on cell and move glial cells appropriately******/
		double forceFact,invPdist2, pdistX,pdistY,pdist2, pdist1;
		int counter1,counter2;	
		for (counter1 = 0; counter1 < numCells -1; counter1++) {
			for (counter2 = counter1+1; counter2 < numCells; counter2++) {
				//Calculate particle-particle distance
				pdistX = cellsX[counter1] - cellsX[counter2];
				pdistY = cellsY[counter1] - cellsY[counter2];

				//Calculate distance squared
				pdist2 = pdistX*pdistX + pdistY*pdistY;

				//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
				//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
				//for details.
				invPdist2 = 1/pdist2;
				forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

				//Calculate the action and reaction for the two particles				
				forcesX[counter1] = forcesX[counter1] - pdistX * forceFact;
				forcesY[counter1] = forcesY[counter1] - pdistY * forceFact;
				forcesX[counter2] = forcesX[counter2] + pdistX * forceFact;
				forcesY[counter2] = forcesY[counter2] + pdistY * forceFact;
			}
			for (counter2 = 0; counter2 < numVesselcells; counter2++) {
				//Calculate particle-particle distance
				pdistX = cellsX[counter1] - VesselcellsX[counter2];
				pdistY = cellsY[counter1] - VesselcellsY[counter2];

				//Calculate distance squared
				pdist2 = pdistX*pdistX + pdistY*pdistY;

				//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
				//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
				//for details.
				invPdist2 = 1/pdist2;
				forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

				//Calculate the action and reaction for the two particles
				forcesX[counter1] = forcesX[counter1] - 2*pdistX * forceFact;
				forcesY[counter1] = forcesY[counter1] - 2*pdistY * forceFact;
			}
		}

		//Update coordinates
		for (counter1 = 0; counter1 < numCells; counter1++) {
			if (forcesX[counter1] > .2) forcesX[counter1] = .2;
			else if (forcesX[counter1] < -.2) forcesX[counter1] = -.2;
			if (forcesY[counter1] > .2) forcesY[counter1] = .2;
			else if (forcesY[counter1] < -.2) forcesY[counter1] = -.2;
			cellsX[counter1] = cellsX[counter1] - 0.5*dt2* 48 *forcesX[counter1];
			cellsY[counter1] = cellsY[counter1]  - 0.5*dt2* 48 *forcesY[counter1];
			//Check if the forces put the cell inside the vessel
			if (cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1] <= (Rvessel+RvesselCells)*(Rvessel+RvesselCells)) {
				if (forcesX[counter1] < 0) cellsX[counter1] = cellsX[counter1] - 2;
				else cellsX[counter1] = cellsX[counter1] +2; 
				if (forcesY[counter1] < 0) cellsY[counter1] = cellsY[counter1] - 2;
				else cellsY[counter1] = cellsY[counter1] +2; 
				cellsX[counter1] = cellsX[counter1] + dt2* 48 *forcesX[counter1];
				cellsY[counter1] = cellsY[counter1] + dt2* 48 *forcesY[counter1];	
			}
			forcesX[counter1] = 0;	
			forcesY[counter1] = 0;	
			cellDist[counter1] = sqrt(cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1])-Rvessel;
		}

		int day = 0; //day of the week in simulation
		int week = 0; //week in simulation
		int hour = 0; //hour in simulation
		int minute = 0; //minute in simulation

		for (T=0; T < therapyStart; T++) {

			/*****Handle Cell Death******/
			int numKillCells = 0;
			int * cellsToKill;
			double t1,t2;
			cellsToKill = (int *) malloc(sizeof(int)*50000);	
			for (i=0; i< numCells; i++)
			{
				cellAge[i] = cellAge[i] +1;
				if (cellDeathTime[i] == 0)
				{
					cellsToKill[numKillCells] = i;
					numKillCells = numKillCells + 1;
					TOT_KILL = TOT_KILL + 1;
				}

			}
			stemCount[T] = 0;
			tbCount[T] = 0;
			for (i=0; i < numCells; i++) {
				if (cellType[i] == 0) {
					stemCount[T] = stemCount[T] + 1;
				}
				else {
					tbCount[T] = tbCount[T] + 1;
				}
			}
			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=cellsToKill[i]; j<numCells-i-1; j++) {				
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellDivTime[numCells] = -1;
					cellConc[numCells] = 0;
					cellDiffTime[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDist[numCells] = 0;
					cellAge[numCells] = -1;
					cellDivRate[numCells] = 0;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}
			free(cellsToKill);
			/*****Handle Cell Division******/
			//Given a Poisson distribution function determining the likelihood of cell division, calculate cell division
			for (i=0; i < numCells; i++) {
				if (cellRadii[i] < R*.1) cellDivTime[i] = (int) floor( randZerotoOne() * (lambda));
				if (T<600) cellDivTime[i] = (int) floor( randZerotoOne() * 180);
				if (cellDivTime[i] == 0)  {
					if (cellQ[i] == 0) {
						TOT_DIV = TOT_DIV + 1;
						cellRadii[numCells] = 0.5*cellRadii[i];
						cellRadii[i] = cellRadii[i] *0.5;
						double rtmp = randZerotoOne();
						if (rtmp < .25) { //put horizontal 	
							cellsX[numCells] = cellsX[i] + cellRadii[i];
							cellsY[numCells] = cellsY[i];
							cellsX[i] = cellsX[i] - cellRadii[i];
						}
						else if (rtmp < 0.5) { //put vertical
							cellsY[numCells] = cellsY[i] + cellRadii[i];
							cellsX[numCells] = cellsX[i];
							cellsY[i] = cellsY[i] - cellRadii[i];
						}
						else if (rtmp < 0.75) { //diagonal
							double h1 = sqrt(cellRadii[i] * cellRadii[i]*0.5);
							cellsX[numCells] = cellsX[i] + h1;
							cellsX[i] = cellsX[i] -h1 ;
							cellsY[numCells] = cellsY[i]+h1;
							cellsY[i] = cellsY[i] -h1;
						}
						else {
							double h1 = sqrt(cellRadii[i] * cellRadii[i]*0.5);
							cellsX[numCells] = cellsX[i] - h1;
							cellsX[i] = cellsX[i] + h1 ;
							cellsY[numCells] = cellsY[i]-h1;
							cellsY[i] = cellsY[i] + h1;
						}

						if (cellType[i] == 0) cellDivRate[i] = (1/rS)*DT*60;
						else cellDivRate[i] = (1/rTB)*DT*60;
						cellDivRate[numCells] = cellDivRate[i];
						cellConc[numCells] = cellConc[i];
						cellDivTime[i] = (int) floor( randZerotoOne() * lambda);
						cellDivTime[numCells] = (int) floor( randZerotoOne() * lambda);
						cellDeathTime[i] = (int) floor( randZerotoOne() * (randomDeathProb));
						cellAge[i]  = cellDeathTime[i];
						cellDeathTime[numCells] = (int) floor( randZerotoOne() * (randomDeathProb));;
						cellAge[numCells]  = cellDeathTime[numCells];
						cellDiffTime[i] = (int) floor( randZerotoOne() * (randomDiff));
						cellDiffTime[numCells] = (int) floor( randZerotoOne() * (randomDiff));;
						celldeDiffTime[i] =  (int) floor( randZerotoOne() * (randomDediff));
						celldeDiffTime[numCells] = -1;
						cellType[numCells] = cellType[i];
						cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])- Rvessel;
						if (cellDist[i] >= 3) {
							if (randZerotoOne() <= .7) cellType[i] = 0;
						}
						else {
							if (randZerotoOne() <= .7) cellType[i] = 1;
						}
						cellDist[numCells] = sqrt(cellsX[numCells]*cellsX[numCells] + cellsY[numCells] * cellsY[numCells]);
						cellDist[numCells] = cellDist[numCells] - Rvessel;
						if (cellDist[i] >= 3) {
							if (randZerotoOne() <= .7) cellType[numCells] = 0;
						}
						else {
							if (randZerotoOne() <= .7) cellType[numCells] = 1;
						}
						cellQ[i] = 0;
						cellQ[numCells] = 0;
						radQ[i] = 0;
						radQt[i] = 0;
						radQ[numCells] = 0;
						radQt[numCells] = 0;
						numCells = numCells + 1;
					}
				}
			}

			/*****Handle Cell Growth******/
			for (i=0; i < numCells; i++) {
				if (cellRadii[i] <= R*1.25) {
					cellRadii[i] = cellRadii[i] * growth;
				}
			}

			/*****Calculate Force on cell and move glial cells appropriately******/
			double forceFact,invPdist2, pdistX,pdistY,pdist2, pdist1;

			for (counter1 = 0; counter1 < numCells -1; counter1++) {
				for (counter2 = counter1+1; counter2 < numCells; counter2++) {
					//Calculate particle-particle distance
					pdistX = cellsX[counter1] - cellsX[counter2];
					pdistY = cellsY[counter1] - cellsY[counter2];

					//Calculate distance squared
					pdist2 = pdistX*pdistX + pdistY*pdistY;

					//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
					//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
					//for details.
					invPdist2 = 1/pdist2;
					forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

					//Calculate the action and reaction for the two particles
					forcesX[counter1] = forcesX[counter1] + pdistX * forceFact;
					forcesY[counter1] = forcesY[counter1] + pdistY * forceFact;
					forcesX[counter2] = forcesX[counter2] - pdistX * forceFact;
					forcesY[counter2] = forcesY[counter2] - pdistY * forceFact;
				}
				for (counter2 = 0; counter2 < numVesselcells; counter2++) {
					//Calculate particle-particle distance
					pdistX = cellsX[counter1] - VesselcellsX[counter2];
					pdistY = cellsY[counter1] - VesselcellsY[counter2];

					//Calculate distance squared
					pdist2 = pdistX*pdistX + pdistY*pdistY;

					//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
					//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
					//for details.
					invPdist2 = 1/pdist2;
					forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

					//Calculate the action and reaction for the two particles
					forcesX[counter1] = forcesX[counter1] + 2*pdistX * forceFact;
					forcesY[counter1] = forcesY[counter1] + 2*pdistY * forceFact;
				}
			}

			//Update coordinates
			for (counter1 = 0; counter1 < numCells; counter1++) {
				if (forcesX[counter1] > .2) forcesX[counter1] = .2;
				else if (forcesX[counter1] < -.2) forcesX[counter1] = -.2;
				if (forcesY[counter1] > .2) forcesY[counter1] = .2;
				else if (forcesY[counter1] < -.2) forcesY[counter1] = -.2;
				cellsX[counter1] = cellsX[counter1] + 0.5*dt2* 48 *forcesX[counter1];
				cellsY[counter1] = cellsY[counter1]  + 0.5*dt2* 48 *forcesY[counter1];
				if (cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1] <= (Rvessel+RvesselCells)*(Rvessel+RvesselCells)) {
					if (forcesX[counter1] < 0) cellsX[counter1] = cellsX[counter1] -2;
					else cellsX[counter1] = cellsX[counter1] +2; 
					if (forcesY[counter1] < 0) cellsY[counter1] = cellsY[counter1] -2;
					else cellsY[counter1] = cellsY[counter1] +2; 
					cellsX[counter1] = cellsX[counter1] + 2*dt2* 48 *forcesX[counter1];
					cellsY[counter1] = cellsY[counter1] + 2*dt2* 48 *forcesY[counter1];	
				}
				forcesX[counter1] = 0;
				forcesY[counter1] = 0;
				cellDist[counter1] = sqrt(cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1])- Rvessel;
			}

			/*****Handle Cell Differentiation******/
			cellsToKill = (int *) malloc(sizeof(int)*50000);
			numKillCells = 0;
			for (i=0; i < numCells; i++) {
				if (cellDiffTime[i] == 0) {
					if (cellQ[i] == 0) {
						cellType[i] = 1; 
						cellDiffTime[i] =  (int) floor(randZerotoOne() * (randomDiff)); 
						cellDivRate[i] = (1.0/nuTumorBulk * DT * 60);
						cellDivTime[i] = mTumorBulk + (int) floor( randZerotoOne() * cellDivRate[i]);
						if (cellType[i] >= Z) {
							cellsToKill[numKillCells] = i;
							numKillCells = numKillCells + 1;
						}
					}
				}
			}
			TOT_KILL = TOT_KILL + numKillCells;

			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=1; j<numCells-i; j++) {
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellDivTime[numCells] = -1;
					cellConc[numCells] = 0;
					cellDiffTime[numCells] = -1;
					cellAge[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDist[numCells] = 0;
					cellDivRate[numCells] = 0;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}

			free(cellsToKill);

			/*****Handle Cell deDifferentiation******/
			for (i=0; i < numCells; i++) {
				if (celldeDiffTime[i] == 0) {
					if (cellQ[i] == 0) {
						cellType[i] = 0;
						celldeDiffTime[i] = (int) floor(randZerotoOne() * (randomDediff));
					}
				}
			}

			/*****Handle Aging******/
			for (i=0; i < numCells; i++) {
				if (cellQ[i] == 0)  {
					cellDivTime[i] = cellDivTime[i] - 1;
					cellDeathTime[i] = cellDeathTime[i] - 1;
					cellDiffTime[i] = cellDiffTime[i] - 1;
					celldeDiffTime[i] = celldeDiffTime[i] - 1;
				}
			}


		}

		for (i=0; i < numCells; i++) {
			if (cellDist[i] < (Rvessel + RvesselCells + tau * R)) cellType[i] = 1;
		}

		stemCount[T] = 0;
		tbCount[T] = 0;
		for (i=0; i < numCells; i++) {
			if (cellType[i] == 0) {
				stemCount[T] = stemCount[T] + 1;
			}
			else {
				tbCount[T] = tbCount[T] + 1;
			}
		}

		int counter = tbCount[T];
		int count2 = 0;
		int tTB = numCells - ceil( (double) numCells / ratioTbtoS );
		for (i=0; i< numCells; i++) {
			if ((counter < tTB)  && (cellDist[i] > 8))
			{
				if (cellType[i] == 0) counter++;  
				cellType[i] = 0;
			}
		}
		stemCount[T] = 0;	
		tbCount[T] = 0;	

		for (i=0; i < numCells; i++) {
			if (cellType[i] == 0) {
				stemCount[T] = stemCount[T] + 1;
			}
			else {
				tbCount[T] = tbCount[T] + 1;
			}
		}

		int ave_death, ave_div;
		ave_death = 0;
		ave_div = 0;
		for (i=0; i< numCells; i++) {
			cellDeathTime[i] = (int) floor(randZerotoOne() * (randomDeathProb));	
			cellDivTime[i] = (int) floor(randZerotoOne() * (lambda));	
			cellAge[i] = cellDeathTime[i];
			ave_death += cellDeathTime[i];
			ave_div += cellDivTime[i];
		}


		for (T=1; T < maxTimeSteps; T++) {
			stemCount[T] = 0;
			tbCount[T] = 0;
			for (i=0; i < numCells; i++) {
				if (cellType[i] == 0) {
					stemCount[T] = stemCount[T] + 1;
				}
				else {
					tbCount[T] = tbCount[T] + 1;
				}
			}

			//Keep Track of what week, day, time it is after therapy commences
			minute = minute + 1;
			if (minute > 60*DT) {
				hour = hour + 1;
				minute = 0;
			}
			if (hour > 24) {
				day = day + 1;
				hour = 0;
			}
			if (day > 7) {
				week = week + 1;
				day = 0;
			}

			/*****Handle Cell Death******/
			int numKillCells = 0;
			int * cellsToKill;
			double t1,t2;
			cellsToKill = (int *) malloc(sizeof(int)*50000);	
			for (i=0; i< numCells; i++)
			{
				if (cellDeathTime[i] == 0)
				{
					cellsToKill[numKillCells] = i;
					numKillCells = numKillCells + 1;
					TOT_KILL = TOT_KILL + 1;
				}

			}
			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=cellsToKill[i]; j<numCells-i-1; j++) {				
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellConc[numCells] = 0;
					cellDiffTime[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDivRate[numCells] = 0;
					cellDist[numCells] = 0;
					cellAge[numCells] = -1;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}
			free(cellsToKill);

			/*****Handle Cell Division******/
			for (i=0; i < numCells; i++) {
				if (cellRadii[i] < R*.1) cellDivTime[i] = (int) floor( randZerotoOne() * (lambda));
				if (cellDivTime[i] == 0)  {
					if (cellQ[i] == 0) {
						TOT_DIV = TOT_DIV + 1;
						cellRadii[numCells] = 0.5*cellRadii[i];
						cellRadii[i] = cellRadii[i] *0.5;
						double rtmp = randZerotoOne();
						if (rtmp < .25) { //put horizontal 	
							cellsX[numCells] = cellsX[i] + cellRadii[i];
							cellsY[numCells] = cellsY[i];
							cellsX[i] = cellsX[i] - cellRadii[i];
						}
						else if (rtmp < 0.5) { //put vertical
							cellsY[numCells] = cellsY[i] + cellRadii[i];
							cellsX[numCells] = cellsX[i];
							cellsY[i] = cellsY[i] - cellRadii[i];
						}
						else if (rtmp < 0.75) { //diagonal
							double h1 = sqrt(cellRadii[i] * cellRadii[i]*0.5);
							cellsX[numCells] = cellsX[i] + h1;
							cellsX[i] = cellsX[i] -h1 ;
							cellsY[numCells] = cellsY[i]+h1;
							cellsY[i] = cellsY[i] -h1;
						}
						else {
							double h1 = sqrt(cellRadii[i] * cellRadii[i]*0.5);
							cellsX[numCells] = cellsX[i] - h1;
							cellsX[i] = cellsX[i] + h1 ;
							cellsY[numCells] = cellsY[i]-h1;
							cellsY[i] = cellsY[i] + h1;
						}
						if (cellType[i] == 0) cellDivRate[i] = (1/rS)*DT*60;
						else cellDivRate[i] = (1/rTB)*DT*60;
						cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])- Rvessel;
						if (cellDist[i] >= 3) {
							if (randZerotoOne() <= .7) cellType[i] = 0;
						}
						else {
							if (randZerotoOne() <= .7) cellType[i] = 1;
						}
						cellDist[numCells] = sqrt(cellsX[numCells]*cellsX[numCells] + cellsY[numCells] * cellsY[numCells]);
						cellDist[numCells] = cellDist[numCells] - Rvessel;
						if (cellDist[i] >= 3) {
							if (randZerotoOne() <= .7) cellType[numCells] = 0;
						}
						else {
							if (randZerotoOne() <= .7) cellType[numCells] = 1;
						}
						cellDivRate[numCells] = cellDivRate[i];
						cellDivTime[i] = (int) floor( randZerotoOne() * lambda);
						cellConc[numCells] = cellConc[i];
						cellDivTime[numCells] = (int) floor( randZerotoOne() * lambda);
						cellDeathTime[i] = (int) floor( randZerotoOne() * (randomDeathProb));
						cellDeathTime[numCells] = (int) floor( randZerotoOne() * (randomDeathProb));
						cellAge[i] = cellDeathTime[i];
						cellAge[numCells] = cellDeathTime[numCells];
						cellDiffTime[i] = (int) floor( randZerotoOne() * (randomDiff));
						cellDiffTime[numCells] = (int) floor( randZerotoOne() * (randomDiff));;
						celldeDiffTime[i] = (int) floor( randZerotoOne() * (randomDediff));
						celldeDiffTime[numCells] = (int) floor( randZerotoOne() * (randomDediff));;
						cellType[numCells] = cellType[i];
						cellQ[i] = 0;
						cellQ[numCells] = 0;
						radQ[i] = 0;
						radQ[numCells] = 0;
						radQt[i] = 0;
						radQt[numCells] = 0;
						numCells = numCells + 1;
					}
				}
			}

			for (i=0; i < numCells; i++) {
				if (cellRadii[i] <= R*1.25) {
					cellRadii[i] = cellRadii[i] * growth;
				}
			}

			/*****Calculate Force on cell and move glial cells appropriately******/
			double forceFact,invPdist2, pdistX,pdistY,pdist2, pdist1;

			for (counter1 = 0; counter1 < numCells -1; counter1++) {
				for (counter2 = counter1+1; counter2 < numCells; counter2++) {
					//Calculate particle-particle distance
					pdistX = cellsX[counter1] - cellsX[counter2];
					pdistY = cellsY[counter1] - cellsY[counter2];

					//Calculate distance squared
					pdist2 = pdistX*pdistX + pdistY*pdistY;

					//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
					//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
					//for details.
					invPdist2 = 1/pdist2;
					forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

					//Calculate the action and reaction for the two particles
					forcesX[counter1] = forcesX[counter1] + pdistX * forceFact;
					forcesY[counter1] = forcesY[counter1] + pdistY * forceFact;
					forcesX[counter2] = forcesX[counter2] - pdistX * forceFact;
					forcesY[counter2] = forcesY[counter2] - pdistY * forceFact;
				}
				for (counter2 = 0; counter2 < numVesselcells; counter2++) {
					//Calculate particle-particle distance
					pdistX = cellsX[counter1] - VesselcellsX[counter2];
					pdistY = cellsY[counter1] - VesselcellsY[counter2];

					//Calculate distance squared
					pdist2 = pdistX*pdistX + pdistY*pdistY;

					//Calculate Lennard-Jones potential assuming sigma=1 and epsilon=1
					//See http://www.cchem.berkeley.edu/chem195/_l_j___force_8m.html#af8855bc03346959adac398ca74c45a06
					//for details.
					invPdist2 = 1/pdist2;
					forceFact = pow(invPdist2,4) * (pow(invPdist2,3)-0.5);

					//Calculate the action and reaction for the two particles
					forcesX[counter1] = forcesX[counter1] + 2*pdistX * forceFact;
					forcesY[counter1] = forcesY[counter1] + 2*pdistY * forceFact;
				}
			}

			//Update coordinates
			for (counter1 = 0; counter1 < numCells; counter1++) {
				if (forcesX[counter1] > .2) forcesX[counter1] = .2;
				else if (forcesX[counter1] < -.2) forcesX[counter1] = -.2;
				if (forcesY[counter1] > .2) forcesY[counter1] = .2;
				else if (forcesY[counter1] < -.2) forcesY[counter1] = -.2;
				cellsX[counter1] = cellsX[counter1] + 0.5*dt2* 48 *forcesX[counter1];
				cellsY[counter1] = cellsY[counter1]  + 0.5*dt2* 48 *forcesY[counter1];
				//Check if the forces put the cell inside the vessel
				if (cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1] <= (Rvessel+RvesselCells)*(Rvessel+RvesselCells)) {
					if (forcesX[counter1] < 0) cellsX[counter1] = cellsX[counter1] -2;
					else cellsX[counter1] = cellsX[counter1] +2; 
					if (forcesY[counter1] < 0) cellsY[counter1] = cellsY[counter1] -2;
					else cellsY[counter1] = cellsY[counter1] +2; 
					cellsX[counter1] = cellsX[counter1] + 2*dt2* 48 *forcesX[counter1];
					cellsY[counter1] = cellsY[counter1] + 2*dt2* 48 *forcesY[counter1];	
				}
				forcesX[counter1] = 0;
				forcesY[counter1] = 0;
				cellDist[counter1] = sqrt(cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1]);
				cellDist[counter1] = cellDist[counter1] - Rvessel;
			}

			cellsToKill = (int *) malloc(sizeof(int)*50000);
			numKillCells = 0;
			for (i=0; i < numCells; i++) {
				if (cellDiffTime[i] == 0) {
					if (cellQ[i] == 0) {
						cellType[i] = cellType[i] + 1;
						cellDiffTime[i] =  (int) floor(randZerotoOne() * (randomDiff)); 
						if (cellType[i] >= Z) {
							cellsToKill[numKillCells] = i;
							numKillCells = numKillCells + 1;
						}
					}
				}
			}
			TOT_KILL = TOT_KILL + numKillCells;

			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=1; j<numCells-i; j++) {
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellDivTime[numCells] = -1;
					cellConc[numCells] =0;
					cellDiffTime[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDist[numCells] = 0;
					cellAge[numCells] = -1;
					cellDivRate[numCells] = 0;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}

			free(cellsToKill);


			/*****Handle Cell deDifferentiation******/
			for (i=0; i < numCells; i++) {
				if (celldeDiffTime[i] == 0) {
					if (cellQ[i] == 0) {
						if (cellType[i] <= Zrevert) {
							cellType[i] = cellType[i] - 1;
							if (cellType[i] < 0 )
								cellType[i] = 0;
							celldeDiffTime[i] = (int) floor(randZerotoOne() * (randomDediff));
						}
					}
				}
			}

			/*****Handle Radiotherapy******/
			giveRad=0; //reset radiation flag
			cellsToKill = (int *) malloc(sizeof(int)*70000);
			numKillCells = 0;
			if (week < RadWeeks) {
				for (i=0; i< radDoses; i++) {
					if (day == radDay[i]) {
						if (hour == radHour[i]) {
							if (minute == 8) { 
								radDose = radDoseA[i];
								if (radDose > 0) { 
									giveRad = 1;
								}
							}
						}
					}
				}
			}
			double tmp5;
			if (giveRad == 1) {
				TOT_RAD = TOT_RAD + 1;
				numKillCells=0;
				for (i=0; i < numCells; i++) {
					if (cellType[i] ==0) {
						rho = RadRho;
						radQ[i] = DT*60*LStem + (int) floor(randZerotoOne() * (1/lambdaStem)*DT*2 );
						radQt[i] = 1;
						alpha = alphaS;
						beta = betaS;
					}
					else {
						rho = 1;
						radQ[i] = DT*60*LTumorBulk + (int) floor(randZerotoOne() * (1/lambdaTumorBulk) *DT*2);
						radQt[i] = 1;
						alpha = alphaT;
						beta = betaT;
					}
					tmp5 = sens*cellConc[i];
					if (tmp5 > 0) tmp5 = tmp5/(.00005);
					else tmp5 = 0;
					if (randZerotoOne() <= 1-exp(-1*alpha*(1+tmp5)*rho*radDose-beta*rho*radDose*radDose)) {
						cellsToKill[numKillCells] = i;
						numKillCells = numKillCells + 1;
						cellType[i] = 3;
						RAD_KILL=RAD_KILL+1;
					}
					else {
						if (timeDepGamma == 1) gamma = .4 * exponential(10,-1*(T-mu)*(T-mu)/(sigmaS));
						celldeDiffTime[i] = -1;
						if (randZerotoOne() <= gamma) {
							if (cellQ[i] == 0) {
								if (cellType[i] < Zrevert) {
									cellType[i] = cellType[i] -1;
									celldeDiffTime[i] = (int) floor(randZerotoOne() * (randomDediff));
									if (cellType[i] < 0) cellType[i] = 0;
								}
							}
						}
					}
				}
			}


			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=cellsToKill[i]; j<numCells-i-1; j++) {				
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellDivTime[numCells] = -1;
					cellConc[numCells] = 0;
					cellDiffTime[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDist[numCells] = 0;
					cellAge[numCells] = -1;
					cellDivRate[numCells] = 0;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}

			for (i=0; i < numCells; i++) {
				if (cellType[i] >1) cellType[i] = 1;
			}
			free(cellsToKill);
			for (i=0; i < numCells; i++) {
				if (radQ[i] > 0) {
					radQ[i] = radQ[i] - 1;
					radQt[i] = radQt[i] + 1;
					cellQ[i] = 1;
					if (radQ[i] ==0) {   
						if (cellType[i] == 0) cellDivRate[i] = (1/rS)*DT*60;
						else cellDivRate[i] = (1/rTB)*DT*60;
						cellDivRate[numCells] = cellDivRate[i];
						cellDivTime[i] = (int) floor( randZerotoOne() * (lambda));
						cellDivTime[numCells] = (int) floor( randZerotoOne() * (lambda));;
					}
				} else {
					radQt[i] = 0;
					cellQ[i] = 0;
				}
			}


			cellsToKill = (int *) malloc(sizeof(int)*50000);
			numKillCells = 0;
			if (ChemoDoses > 0) {
				if (week < ChemoWeeks) {
					for (i=0; i< ChemoDoses; i++) {
						if (day == ChemoDay[i]) {
							if (hour == ChemoHour[i]) {
								if (minute == 0) {
									ChemoTime[chemoDosesGiven] = T;
									ThalfPerDose[chemoDosesGiven] = T + Thalf; 
									chemoDosesGiven = chemoDosesGiven + 1;
								}
							}
						}
					}
				}	
			}
			double r,r2;

			for (i=0; i < chemoDosesGiven; i++) {
				if (T < ChemoTime[i] + Tmax*DT) {
					r = -1*(T-ChemoTime[i]) / (double) Thalf;
					r2 = pow(2,r);
					r = Cmax/r2+1;
					r2=log(r);
					r=r2/(60.0);
					r2 = exp(r*(T-ChemoTime[i]));
					r = pow(2,-1*(T-ChemoTime[i])/ (double) Thalf);
					ChemoC[i] = r*(r2-1);
					if (ChemoC[i] < 0) ChemoC[i] =0;
				} else {
					ChemoC[i] = Cmax * ((pow(2,-1*((T-ChemoTime[i] - Tmax*DT)/(double) Thalf))));
				}
			}

			double F, cCell = 0, tau=0, Ttemp=0;
			if (chemoDosesGiven >0) {
				for (i=0; i < numCells; i++) {
					cCell = 0;
					for (j=0; j < chemoDosesGiven; j++) {
						Ttemp = 2 * Thalf;
						tau = T - ChemoTime[j];
						if ((cellDist[i] < 15) && (T-ChemoTime[j] < 2*Thalf)) cCell = cCell + ChemoC[j]*(1-cellDist[i]/15.0);
						cellConc[i] = cCell;

					}
					F = cellConc[i]/IC50;//100 * cCell * IC50 / (1+IC50*cCell);
					if (randZerotoOne() <= F) {
						cellsToKill[numKillCells] = i;
						numKillCells = numKillCells + 1;
						Ch_rm=Ch_rm+1;
						cellType[i] = 3;
					} 
				}
			}
			if (numKillCells > 0) {
				for (i=0; i<numKillCells; i++) {
					for (j=cellsToKill[i]; j<numCells-i-1; j++) {				
						cellsX[j] = cellsX[1+j];
						cellsY[j] = cellsY[1+j];
						cellDivTime[j] = cellDivTime[1+j];
						cellConc[j] = cellConc[1+j];
						cellDeathTime[j] = cellDeathTime[1+j];
						cellDiffTime[j] = cellDiffTime[1+j];
						celldeDiffTime[j] = celldeDiffTime[1+j];
						cellType[j] = cellType[1+j];
						cellDist[j] = cellDist[1+j];
						cellAge[j] = cellAge[1+j];
						cellDivRate[j] = cellDivRate[1+j];
						cellQ[j] = cellQ[1+j];
						radQ[j] = radQ[1+j];
						radQt[j] = radQt[1+j];
						cellRadii[j] = cellRadii[1+j];
					}
					cellsX[numCells] = 0;
					cellsY[numCells] = 0;
					cellDeathTime[numCells] = -1;
					cellDivTime[numCells] = -1;
					cellConc[numCells] = 0;
					cellDiffTime[numCells] = -1;
					celldeDiffTime[numCells] = -1;
					cellType[numCells] = 0;
					cellDist[numCells] = 0;
					cellAge[numCells] = -1;
					cellDivRate[numCells] = 0;
					cellQ[numCells] = 0;
					cellRadii[numCells] = 0;
					radQ[numCells] = 0;
					radQt[numCells] = 0;
					numCells = numCells -1;
				}
			}

			for (i=0; i < numCells; i++) {
				if (cellType[i] >1) cellType[i] = 1;
			}
			free(cellsToKill);


			/*****Handle Stem Cell Retention Near Vessel Boundary******/
			for (i=0; i < numCells; i++) {
				if (cellDist[i] < (Rvessel + RvesselCells + 3 * R)) cellType[i] = 1;
				else {
					if (cellType[i] == 0) cellType[i] = 0;
				}
			}
			/*****Handle Aging******/
			for (i=0; i < numCells; i++) {
				if (cellQ[i] == 0)  {
					cellDivTime[i] = cellDivTime[i] - 1;
					cellDeathTime[i] = cellDeathTime[i] - 1;
					cellDiffTime[i] = cellDiffTime[i] - 1;
					celldeDiffTime[i] = celldeDiffTime[i] - 1;
				}
			}

			/*****Count Cell Type******/
			stemCount[T] = 0;
			tbCount[T] = 0;
			for (i=0; i < numCells; i++) {
				if (cellType[i] == 0) {
					stemCount[T] = stemCount[T] + 1;
				}
				else {
					tbCount[T] = tbCount[T] + 1;
				}
			}
		}//timesteps
	}//cycle
}
