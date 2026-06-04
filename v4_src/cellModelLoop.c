#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "string.h"
#include "math.h"



static uint64_t rng_state;

static inline void seed_rng(uint64_t seed)
{
    rng_state = seed ? seed : 1ULL;
}

static inline uint64_t splitmix64(void)
{
    uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline double randZerotoOne(void)
{
    return (splitmix64() >> 11) * (1.0 / 9007199254740992.0);
}

static int cmp_int_asc(const void *a, const void *b)
{
	int ia = *(const int*)a;
	int ib = *(const int*)b;
	return (ia > ib) - (ia < ib);
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

static inline double steps_to_hours(int steps, int DT) {
	return (double)steps / (60.0 * (double)DT);
}

static inline double Cblood_TMQ(double t_hr,
		double Cmax,
		double Thalf_hr,
		double Tmax_hr)
{
	if (t_hr <= 0.0) return 0.0;

	const double ln2 = log(2.0);
	const double decay_to_tmax = exp(-ln2 * (Tmax_hr / Thalf_hr)); // 2^(-Tmax/Thalf)

	// r = ln( Cmax / 2^(-Tmax/Thalf) + 1 ) / Tmax
	const double r = log(Cmax / decay_to_tmax + 1.0) / Tmax_hr;

	if (t_hr < Tmax_hr) {
		const double rise = exp(r * t_hr) - 1.0;
		const double decay = exp(-ln2 * (t_hr / Thalf_hr));       // 2^(-t/Thalf)
		double C = rise * decay;
		return (C > 0.0 ? C : 0.0);
	} else {
		const double decay = exp(-ln2 * ((t_hr - Tmax_hr) / Thalf_hr));
		double C = Cmax * decay;
		return (C > 0.0 ? C : 0.0);
	}
}


static inline int exp_wait_steps(double mean_hours, int DT)
{
	if (mean_hours <= 0.0) return 0;
	double U = randZerotoOne();
	if (U < 1e-12) U = 1e-12;
	double dur_hours = -mean_hours * log(U);
	double steps = dur_hours * (60.0 * (double)DT);
	return (int) llround(steps);
}


/* 
 * Generic remove function: 
 * arrays are passed as references (pointers to pointers),
 * so this function just operates on them in place.
 */
void remove_dead_cells(int * hasExitedQ, int *numCells, double *cellsX, double *cellsY, int *cellDeathTime, int *cellDivTime, double *cellConc, int *cellDiffTime, int *celldeDiffTime, int *cellType, double *cellDist, int *cellAge,  double *cellRadii, int *cellQ, int *radQ, int *radQt, int *cellsToKill, int numKillCells)
{
	if (*numCells <= 0 || numKillCells <= 0) return;

	int n = 0;
	for (int k = 0; k < numKillCells; ++k) {
		int idx = cellsToKill[k];
		if (idx >= 0 && idx < *numCells) {
			cellsToKill[n++] = idx;
		}
	}
	if (n <= 0) return;

	qsort(cellsToKill, n, sizeof(int), cmp_int_asc);

	int u = 1;
	for (int k = 1; k < n; ++k) {
		if (cellsToKill[k] != cellsToKill[u - 1]) {
			cellsToKill[u++] = cellsToKill[k];
		}
	}
	for (int k = u - 1; k >= 0; --k) {
		int dead = cellsToKill[k];
		if (dead < 0 || dead >= *numCells) continue;

		for (int j = dead; j < *numCells - 1; ++j) {
			cellsX[j]         = cellsX[j + 1];
			cellsY[j]         = cellsY[j + 1];
			cellDeathTime[j]  = cellDeathTime[j + 1];
			cellDivTime[j]    = cellDivTime[j + 1];
			hasExitedQ[j]     = hasExitedQ[j + 1];
			cellConc[j]       = cellConc[j + 1];
			cellDiffTime[j]   = cellDiffTime[j + 1];
			celldeDiffTime[j] = celldeDiffTime[j + 1];
			cellType[j]       = cellType[j + 1];
			cellDist[j]       = cellDist[j + 1];
			cellAge[j]        = cellAge[j + 1];
			cellQ[j]          = cellQ[j + 1];
			radQ[j]           = radQ[j + 1];
			radQt[j]          = radQt[j + 1];
			cellRadii[j]      = cellRadii[j + 1];
		}

		int last = *numCells - 1;
		cellsX[last]         = 0.0;
		cellsY[last]         = 0.0;
		cellDeathTime[last]  = -1;
		cellDivTime[last]    = -1;
		hasExitedQ[last]     = -1;
		cellConc[last]       = 0.0;
		cellDiffTime[last]   = -1;
		celldeDiffTime[last] = -1;
		cellType[last]       = 1;
		cellDist[last]       = 0.0;
		cellAge[last]        = -1;
		cellQ[last]          = 0;
		radQ[last]           = 0;
		radQt[last]          = 0;
		cellRadii[last]      = 0.0;

		(*numCells)--;
	}
}



void cellModelLoop(int run_id, int maxTimeSteps, double sens, int Z, int Zrevert, int ChemoWeeks, int ChemoDays, int ChemoDoses, int ThalfPerDose[7], int ChemoDay[7], int ChemoHour[7], double ChemoC[7], double Cmax, int RadWeeks, int radDoses, int radDose, int radDoseA[10], int radHour[10], int radDay[10], int *stemCount, int * tbCount, int *cellAge)
{
	int * cellsToKill = (int *) malloc(sizeof(int)*50000);	
	int age1Counter, age2Counter,typeCounter;
	double dist1Counter, dist2Counter, maxDist;
	char outputFileName[256];
	char time3outputFileName[256];
	char time1outputFileName[256];
	outputFileName[0] = 0;
	time3outputFileName[0] = 0;
	time1outputFileName[0] = 0;

	sprintf(outputFileName, "output_%d.csv", run_id);
	sprintf(time3outputFileName, "output_all_%d.csv", run_id);
	sprintf(time1outputFileName, "celldata_%d.csv", run_id);

	seed_rng((uint64_t)run_id);

	/******Variable Declarations*******/
	int DT = 2; //time steps per minute
	FILE *fp,*fp1, *fp3;
	int ratioTbtoS = 20; //ratio of tumor bulk to stem cells before therapy
	double tau = 2; //number of cell diameters away from vessel wall that stem cells can exit
	int RvesselCells=2; //radius of the cells making up the blood vessel
	int Rvessel = 4; //radius of the blood vessel
	int R=2; //radius of glial cells
	int k = 5000; //spring constant
	const int sens_on = (sens != 0.0);   // sens is a 0/1 flag
	const double CHEMO_RADIUS_CELLS = 8.0;
	double rTB = .0038; //from Cell paper, proliferation rate of tumor cells after exiting quiescence
	double rS = .008; //from Cell paper, proliferation rate of stem cells after exiting quiescence

	double expected_kills = 0.0;
	int stemN = 0, tbN = 0;
	int timeDepGamma = 1; //flag to use time dependent gamma
	double mu = 3.25*DT*60; //time to perak reverstion after irradiation
	double sigmaS = pow(1.46 * 60.0 * DT, 2.0); // variance in timestep units
	double lambda = 8*24*60*DT; 
	double randomDeathProb=(36*24*60*DT); //Probability of cell randomly dying in any given minute
	double aDiff = .0019; //rate at which stem cells convert to tumor bulk
	double aDeDiff = .45; //rate at which tumor bulk convert to stem cells
	double randomDiff= 1/aDiff *DT *60; //Probability of randomly differentiating and dropping a level in the differentiation cascade
	double randomDediff=1/aDeDiff *DT *60;//Probability of randomly dedifferentiating and gaining a level in the differentiation cascade
	double Diffusion = .00066;//Diffusion constant based on molecular weight of 194.15
				  //Citation: Baer et al, Depletion of 06-alkylguanine-DNA alkyl transferase correlates with potentiation of temozolomide and CCNU toxicity in human tumour cells
	double IC50 = .007668;//unit: mol/m3, micromolar to get the IC50 of 22 micromolar seen after 3 hour drug incubation for U373MG in:
			      // Wedge, Stephen R., et al. "In vitro evaluation of temozolomide combined with X-irradiation." Anti-cancer drugs 8.1 (1997): 92-97.

			      //Treatment Parameters
	double alphaS = .00987; //fit from figure 1D in Cell paper, likelihood of death at a single DNA strand break for stem
	double alphaT = .0987; //fit from figure 1D in Cell paper, likelihood of death at a single DNA strand break for TB
	double betaS = 1.14*.00000001; //fit from figure 1D in Cell paper, likelihood of death at a double DNA strand break for stem
	double betaT = 1.14*.0000001; //fit from figure 1D in Cell paper, likelihood of death at a double DNA strand break for TB
	double gamma = .4; //Cell paper Table 2, probability of dedifferentiation due to radiation leading to side population cells
	double RadRho = .4; //Cell Table 2
	double LStem = 36; //hours of halting for stem cells after radiation
	double LTumorBulk = 24; //hours of halting for tumor bulk cells after  radiation
	double lambdaStem = .00328; //exponential rate at which stem cells exit quiescence
	double lambdaTumorBulk =.1; //exponential rate at which tumor bulk cells quit quiescence
	double mTumorBulk = 366*(DT*60); //Minimum time for newly converted stem cells to begin clonal expansion
	double nuTumorBulk = .054; //Rate at which newly convered DSC lead to clonal expansion

	int chemoDosesGiven = 0; //tracks number of chemo doses that have been administered
	const int steps_per_hour = 60 * DT;

	double Tmax_hours  = 0.50;  // mouse example

	for (int j = 0; j < 7; ++j) {
		ChemoC[j]    = 0.0;     // no drug before dosing
	}

	const double hours_per_step = 1.0 / (double)steps_per_hour;

	// arrays
	int ChemoTime[7];      
	for (int j = 0; j < 7; ++j) ChemoTime[j] = -1;


	int giveRad = 0; //binary flag to determine if radiation is applied in this time step
	int RAD_KILL=0;
	int lastRadT = -60000000;
	int Ch_rm=0;
	int TOT_RAD=0;
	int TOT_DIV = 0;
	int TOT_KILL = 0;

	double ct1= 1; 
	double dt2 = ct1/((double) DT) * ct1/((double) DT);
	int therapyStart = 2800;
	int T = 1; //Time in minutes

	double growth =  pow(2.0, 1.0 / (3.0 * 24.0 * 60.0 * (double)DT)); 
	for (T=0; T < maxTimeSteps; T++) {
		//Initialize cell counts for this time step
		tbCount[T] = 0;
		stemCount[T] = 0;
	}

	double cellRadii[60000];
	double cellsX[60000];
	double cellsY[60000];
	double forcesX[60000];
	double forcesY[60000];
	int cellQ[60000]; //this variable determines if this cell is in quiescence or not
	int cellDivTime[60000]; //this variable measures how far into the cell cycle
	int hasExitedQ[60000]; //this variable tracks if that cell has exited quiescence or not 
	double cellConc[60000]; //this variable tracks chemo concentration at that cell
	int cellDeathTime[60000]; //this variable measures how far into the cell cycle
	int cellDiffTime[60000]; //this variable measures how far into the cell cycle
	int celldeDiffTime[60000]; //this variable measures how far into the cell cycle
	int cellType[60000]; //this determine what type the cell is, 0 is for BSTC, 1 is tumor bulk
	double cellDist[60000]; //this determine how far from the vessel wall the center of the cell is
	int radQ[60000]; //tracking time in radiation enforced quiescence
	int radQt[60000]; //tracking time in radiation enforced quiescence timing
	int i, j;
	int numCells = floor( 3.14*2*(Rvessel+2*R)/(2*R)); //number of glial cells to start with

	int max_cycles=1;
	int cycle;
	for (cycle = 0; cycle < max_cycles; cycle++) {
		for (i=0; i < 60000; i++)
		{	
			cellRadii[i] = R;
			cellsX[i] = 0;
			forcesX[i] = 0;
			cellsY[i] = 0;
			forcesY[i] = 0;
			cellQ[i] = 0;
			cellConc[i] = 0;
			hasExitedQ[i] = 0;
			cellDivTime[i] = -1; //60000*2;
			cellDeathTime[i] = -1; //60000*2;
			celldeDiffTime[i] = -1; //60000*2;
			cellDiffTime[i] = -1; //60000*2;
			cellType[i] = 1;
			cellDist[i] = 0;
			cellAge[i] = 0;
			radQ[i] = 0;
			radQt[i] = 0;
			hasExitedQ[i] = 0;
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
			cellDivTime[i] = (int) floor(randZerotoOne() * 120); 
			hasExitedQ[i]=0;
			cellConc[i] = 0;
			cellDeathTime[i] = (int) floor(randZerotoOne() * (randomDeathProb));
			celldeDiffTime[i] = -1; 
			cellDiffTime[i] = -1;
			cellAge[i] = cellDeathTime[i];
			cellType[i] = 0;
			cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])-Rvessel;
			if (cellType[i] == 0) stemCount[0] = stemCount[0] + 1;
			else tbCount[0]++;
			if (cellDeathTime[i] < cellDivTime[i]) countDeath ++;
			else countDiv++;
		}


		/*****Calculate Force on cell and move glial cells appropriately******/
		double forceFact, pdistX,pdistY,pdist2, pdist1;
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
				const double eps = 1e-12;
				if (pdist2 < eps) pdist2 = eps;
				double inv2 = 1.0 / pdist2;      
				double inv4 = inv2 * inv2;       
				double inv6 = inv4 * inv2;       
				double inv8 = inv4 * inv4;       
				forceFact = inv8 * (inv6 - 0.5);

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
				const double eps = 1e-12;
				if (pdist2 < eps) pdist2 = eps;
				double inv2 = 1.0 / pdist2;      
				double inv4 = inv2 * inv2;      
				double inv6 = inv4 * inv2;       
				double inv8 = inv4 * inv4;       
				forceFact = inv8 * (inv6 - 0.5);

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
		int t_minutes = T / DT;                // integer minutes since therapy start
		int minute_of_hour = t_minutes % 60;
		hour = (t_minutes / 60) % 24;
		day  = (t_minutes / (60*24)) % 7;
		week =  t_minutes / (60*24*7);


		fp = fopen("std_output.txt", "w");
		for (T=0; T < therapyStart; T++) {
			int t_minutes = T / DT;
			int minute_of_hour = t_minutes % 60;
			int hour = (t_minutes / 60) % 24;
			int day  = (t_minutes / (60*24)) % 7;
			int week =  t_minutes / (60*24*7);
			/*****Handle Cell Death******/
			int numKillCells = 0;
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
				remove_dead_cells(hasExitedQ, &numCells, cellsX, cellsY, cellDeathTime, cellDivTime, cellConc, cellDiffTime, celldeDiffTime, cellType, cellDist, cellAge, cellRadii, cellQ, radQ, radQt, cellsToKill, numKillCells);

			}

			/*****Handle Cell Division******/
			int cellCount = numCells;
			//Given a Poisson distribution function determining the likelihood of cell division, calculate cell division
			for (i=0; i < cellCount; i++) {
				if (cellRadii[i] < R*.1) {
					cellDivTime[i] =  (int)floor(randZerotoOne() * lambda);
				}
				// in pre-therapy division logic
				if (T < 600 ) {
					cellDivTime[i] = exp_wait_steps(18.0, DT);  // mean 18 hr
				}
				if (cellDivTime[i] <= 0)  {
					if (cellQ[i] <= 0) {
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

						cellConc[numCells] = cellConc[i];
						double rate_i  = (cellType[i] == 0) ? rS : rTB;                 // /hour
						double mean_i  = (hasExitedQ[i] ? ((rate_i > 0.0) ? (1.0/rate_i) : 0.0)
								: 96.0);
						cellDivTime[i] = exp_wait_steps(mean_i, DT);
						int minDivSteps = (int) llround(((cellType[i] == 0) ? 8.0 : 4.0) * 60.0 * DT);
						cellDivTime[i] += minDivSteps;
						minDivSteps = (int) llround(((cellType[numCells] == 0) ? 8.0 : 4.0) * 60.0 * DT);
						cellDivTime[numCells] += minDivSteps;

						double rate_n  = (cellType[numCells] == 0) ? rS : rTB;
						double mean_n  = (hasExitedQ[numCells] ? ((rate_n > 0.0) ? (1.0/rate_n) : 0.0)
                                                                : 96.0);
						cellDivTime[numCells] = exp_wait_steps(mean_n, DT);

						cellDeathTime[i] = (int) floor( randZerotoOne() * (randomDeathProb));
						cellAge[i]  = cellDeathTime[i];
						cellDeathTime[numCells] = (int) floor( randZerotoOne() * (randomDeathProb));;
						cellAge[numCells]  = cellDeathTime[numCells];
						celldeDiffTime[i] =  -1;
						celldeDiffTime[numCells] = -1;
						cellDiffTime[i] = -1;  // parent GSC stays stem-like
						cellDiffTime[numCells] = -1;
						hasExitedQ[numCells] = 1;
						cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])- Rvessel;
						cellDist[numCells] = sqrt(cellsX[numCells]*cellsX[numCells] + cellsY[numCells] * cellsY[numCells]);
						cellDist[numCells] = cellDist[numCells] - Rvessel;
						cellQ[i] = 0;
						cellQ[numCells] = 0;
						radQ[i] = 0;
						radQt[i] = 0;
						radQ[numCells] = 0;
						radQt[numCells] = 0;

						if (cellType[i] == 0) {
							// GSC division: parent remains GSC.
							// Daughter is usually DSC/DTC; rarely self-renews as GSC.
							const double p_asym_division = 0.98;

							cellType[i] = 0;

							if (randZerotoOne() < p_asym_division) {
								cellType[numCells] = 1;
								// Newly created DSC undergoes a maturation delay before clonal expansion
								cellDiffTime[numCells] =
									mTumorBulk + exp_wait_steps(1.0 / nuTumorBulk, DT);
							} else {
								cellType[numCells] = 0;
								cellDiffTime[numCells] = -1;
							}

						} else {
							// DTC division: both remain DTC
							cellType[i] = 1;
							cellType[numCells] = 1;
						}
						numCells = numCells + 1;
						if (numCells > 59900) {
							printf("TOO MANY CELLS! %d %d %d \n", T, stemCount[T-1], tbCount[T-1]);
							abort();	
						}
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
			double forceFact, pdistX,pdistY,pdist2, pdist1;

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
					const double eps = 1e-12;
					if (pdist2 < eps) pdist2 = eps;
					double inv2 = 1.0 / pdist2;      
					double inv4 = inv2 * inv2;       
					double inv6 = inv4 * inv2;       
					double inv8 = inv4 * inv4;       
					forceFact = inv8 * (inv6 - 0.5);

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
					const double eps = 1e-12;
					if (pdist2 < eps) pdist2 = eps;
					double inv2 = 1.0 / pdist2;      
					double inv4 = inv2 * inv2;       
					double inv6 = inv4 * inv2;       
					double inv8 = inv4 * inv4;       
					forceFact = inv8 * (inv6 - 0.5);

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
				cellDist[counter1] = sqrt(cellsX[counter1]*cellsX[counter1] + cellsY[counter1] * cellsY[counter1])- Rvessel;
			}

			/*****Handle Cell Differentiation******/
			numKillCells = 0;
			for (i = 0; i < numCells; i++) {
				if (cellDiffTime[i] == 0) {
					if (cellQ[i] <= 0) {

						cellType[i] = 1; // becomes tumor bulk
						cellDiffTime[i] = -1;


						if (hasExitedQ[i] == 0) {
							// Baseline division: mean 24 hours
							double rate = (cellType[i] == 0) ? rS : rTB;
							double mean = (hasExitedQ[i] ? ((rate > 0.0) ? (1.0/rate) : 0.0)
									:96.0);
							cellDivTime[i] = mTumorBulk + exp_wait_steps(mean, DT);
						} else {
							// Post-quiescence proliferative division: rS/rTB are rates (/hour) -> convert to mean hours
							double rate = (cellType[i] == 0) ? rS : rTB;
							double mean = (hasExitedQ[i] ? ((rate > 0.0) ? (1.0/rate) : 0.0)
									:96.0);
							cellDivTime[i] = mTumorBulk + exp_wait_steps(mean, DT);
						}

						if (cellType[i] >= Z) {
							cellsToKill[numKillCells] = i;
							numKillCells++;
						}
					}
				}
			}

			TOT_KILL = TOT_KILL + numKillCells;

			if (numKillCells > 0) {
				remove_dead_cells(hasExitedQ, &numCells, cellsX, cellsY, cellDeathTime, cellDivTime, cellConc, cellDiffTime, celldeDiffTime, cellType, cellDist, cellAge, cellRadii, cellQ, radQ, radQt, cellsToKill, numKillCells);
			}


			/*****Handle Cell deDifferentiation******/
			/* Only cells explicitly marked for reversion are allowed to dedifferentiate.
			 * Radiation block should set celldeDiffTime[i] >= 0 only for eligible DTCs.
			 */
			for (i = 0; i < numCells; i++) {
				if (cellQ[i] <= 0 &&
						cellType[i] == 1 &&
						celldeDiffTime[i] == 0) {

					cellType[i] = 0;  // DTC -> GSC

					// Turn reversion timer off until another radiation-triggered event
					celldeDiffTime[i] = -1;
				}
			}
			/*****Handle Aging******/
			for (i=0; i < numCells; i++) {
				cellDeathTime[i] = cellDeathTime[i] - 1;

				if (cellQ[i] == 0)  {
					cellDivTime[i] = cellDivTime[i] - 1;
					if (cellDiffTime[i] > 0) cellDiffTime[i] = cellDiffTime[i] - 1;
					if (celldeDiffTime[i] > 0) celldeDiffTime[i] = celldeDiffTime[i] - 1;
				}

				if (cellDivTime[i] < 0) cellDivTime[i] = 0;
			}


		}

		/*****Handle Stem Cell Retention Near Vessel Boundary******/
		for (i=0; i < numCells; i++) {
			if (cellDist[i] < (tau * R)) {
				if (cellType[i] == 1 && randZerotoOne() < 0.10) {
					cellType[i] = 0;
				}
			}		
		}

		/*****Count Cell Type******/
		stemCount[0] = 0;
		tbCount[0] = 0;
		for (i=0; i < numCells; i++) {
			if (cellType[i] == 0) {
				stemCount[0] = stemCount[0] + 1;
			}
			else {
				tbCount[0] = tbCount[0] + 1;
			}
		}

		/*****Handle Stem Cell to TB Ratio Enforcement******/
		int counter = tbCount[0];
		int targetTB = numCells - (numCells+ratioTbtoS)/(ratioTbtoS+1);
		for (i=0; i< numCells; i++) {
			if ((counter < targetTB)  && (cellDist[i] > 8))
			{
				if (cellType[i] == 0) counter++;  //check if it isn't already set to 1	
				cellType[i] = 1;
			}
		}
		stemCount[0] = 0;	
		tbCount[0] = 0;	

		/*****Count Cell Type******/
		for (i=0; i < numCells; i++) {
			if (cellType[i] == 0) {
				stemCount[0] = stemCount[0] + 1;
			}
			else {
				tbCount[0] = tbCount[0] + 1;
			}
		}


		//Reset div and death time to be even and stable going into treatment
		for (i = 0; i < numCells; i++) {
			cellDeathTime[i] = (int) floor(randZerotoOne() * randomDeathProb);

			double rate = (cellType[i] == 0) ? rS : rTB;  // /hour
			double mean = hasExitedQ[i] ? ((rate > 0.0) ? (1.0 / rate) : 0.0)
				: 96.0;

			cellDivTime[i] = exp_wait_steps(mean, DT);   
			cellAge[i] = cellDeathTime[i];
		}

		fp3 = fopen(time3outputFileName, "w");
		fp1 = fopen(time1outputFileName, "w");

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
			t_minutes = T/DT;
			minute_of_hour = t_minutes % 60;
			hour = (t_minutes / 60) % 24;
			day  = (t_minutes / (60*24)) % 7;
			week =  t_minutes / (60*24*7);

			/*****Handle Cell Death******/
			int numKillCells = 0;
			for (i=0; i< numCells; i++)
			{
				if (cellDeathTime[i] <= 0)
				{
					cellsToKill[numKillCells] = i;
					numKillCells = numKillCells + 1;
					TOT_KILL = TOT_KILL + 1;
				}

			}
			if (numKillCells > 0) {
				remove_dead_cells(hasExitedQ, &numCells, cellsX, cellsY, cellDeathTime, cellDivTime, cellConc, cellDiffTime, celldeDiffTime, cellType, cellDist, cellAge, cellRadii, cellQ, radQ, radQt, cellsToKill, numKillCells);
			}

			/*****Handle Cell Division******/
			int cellCount = numCells;
			//Given a Poisson distribution function determining the likelihood of cell division, calculate cell division
			for (i=0; i < cellCount; i++) {
				if (cellRadii[i] < R*.1) cellDivTime[i] = (int) floor( randZerotoOne() * (lambda));
				if (cellDivTime[i] <= 0)  {
					if (cellQ[i] <= 0) {
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
						cellDist[i] = sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])- Rvessel;
						cellDist[numCells] = sqrt(cellsX[numCells]*cellsX[numCells] + cellsY[numCells] * cellsY[numCells]);
						cellDist[numCells] = cellDist[numCells] - Rvessel;

						double rate = (cellType[i] == 0) ? rS : rTB;
						double mean = (hasExitedQ[i] ? ((rate > 0.0) ? (1.0 / rate) : 0.0) : 96.0);
						cellDivTime[i] = exp_wait_steps(mean, DT);

						rate = (cellType[numCells] == 0) ? rS : rTB;
						mean = (hasExitedQ[numCells] ? ((rate > 0.0) ? (1.0 / rate) : 0.0) : 96.0);
						cellDivTime[numCells] = exp_wait_steps(mean, DT);

						int minDivSteps = (int) llround(((cellType[i] == 0) ? 8.0 : 4.0) * 60.0 * DT);
						cellDivTime[i] += minDivSteps;
						minDivSteps = (int) llround(((cellType[numCells] == 0) ? 8.0 : 4.0) * 60.0 * DT);
						cellDivTime[numCells] += minDivSteps;
						cellConc[numCells] = cellConc[i];
						cellDeathTime[i] = (int) floor( randZerotoOne() * (randomDeathProb));
						cellDeathTime[numCells] = (int) floor( randZerotoOne() * (randomDeathProb));
						hasExitedQ[numCells] = 1;
						cellAge[i] = cellDeathTime[i];
						cellAge[numCells] = cellDeathTime[numCells];
						cellDiffTime[i] = -1; 
						cellDiffTime[numCells] = -1;
						celldeDiffTime[i] = -1;
						celldeDiffTime[numCells] = -1;

						if (cellType[i] == 0) {
							// GSC division: parent remains GSC.
							// Daughter is usually DSC/DTC; rarely self-renews as GSC.
							const double p_asym_division = 0.98;

							cellType[i] = 0;
							cellDiffTime[i] = -1;

							if (randZerotoOne() < p_asym_division) {
								cellType[numCells] = 1;

								// Newly created DSC has delayed clonal expansion
								cellDiffTime[numCells] =
									mTumorBulk + exp_wait_steps(1.0 / nuTumorBulk, DT);
							} else {
								cellType[numCells] = 0;
								cellDiffTime[numCells] = -1;
							}

						} else {
							// DTC division: both remain DTC
							cellType[i] = 1;
							cellType[numCells] = 1;

							cellDiffTime[i] = -1;
							cellDiffTime[numCells] = -1;
						}

						cellQ[i] = 0;
						cellQ[numCells] = 0;
						radQ[i] = 0;
						radQ[numCells] = 0;
						radQt[i] = 0;
						radQt[numCells] = 0;
						numCells = numCells + 1;
						if (numCells > 59900) {
							printf("TOO MANY CELLS! %d %d %d \n", T, stemCount[T-1], tbCount[T-1]);
							abort();	
						}
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
			double forceFact, pdistX,pdistY,pdist2, pdist1;

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
					const double eps = 1e-12;
					if (pdist2 < eps) pdist2 = eps;
					double inv2 = 1.0 / pdist2;      // inv r^2
					double inv4 = inv2 * inv2;       // inv r^4
					double inv6 = inv4 * inv2;       // inv r^6
					double inv8 = inv4 * inv4;       // inv r^8
					forceFact = inv8 * (inv6 - 0.5);

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
					const double eps = 1e-12;
					if (pdist2 < eps) pdist2 = eps;
					double inv2 = 1.0 / pdist2;      
					double inv4 = inv2 * inv2;       
					double inv6 = inv4 * inv2;       
					double inv8 = inv4 * inv4;       
					forceFact = inv8 * (inv6 - 0.5);

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

			/*****Handle Cell Differentiation******/
			numKillCells = 0;
			for (i=0; i < numCells; i++) {
				if (cellDiffTime[i] == 0) {
					if (cellQ[i] == 0) {
						cellType[i] = 1;  //remains tumor bulk after maturation
						cellDiffTime[i] =  -1;
						if (cellType[i] >= Z) {
							//Becomes terminally differentiated and initiate cell death
							cellsToKill[numKillCells] = i;
							numKillCells = numKillCells + 1;
						}
					}
				}
			}
			TOT_KILL = TOT_KILL + numKillCells;

			if (numKillCells > 0) {
				remove_dead_cells(hasExitedQ, &numCells, cellsX, cellsY, cellDeathTime, cellDivTime, cellConc, cellDiffTime, celldeDiffTime, cellType, cellDist, cellAge,  cellRadii, cellQ, radQ, radQt, cellsToKill, numKillCells);
			}



			/*****Handle Cell deDifferentiation******/
			/* Only cells explicitly marked for reversion are allowed to dedifferentiate.
			 * Radiation block should set celldeDiffTime[i] >= 0 only for eligible DTCs.
			 */
			for (i = 0; i < numCells; i++) {
				if (cellQ[i] <= 0 &&
						cellType[i] == 1 &&
						celldeDiffTime[i] == 0) {

					cellType[i] = 0;  // DTC -> GSC

					// Turn reversion timer off until another radiation-triggered event
					celldeDiffTime[i] = -1;
				}
			}

			/*****Handle Chemotherapy****/
			numKillCells = 0;
			//Determine if Dose is administered
			if (ChemoDoses > 0) {
				if (week < ChemoWeeks) {
					for (i=0; i< ChemoDoses; i++) {
						if (day == ChemoDay[i] && hour == ChemoHour[i] && minute_of_hour == 0 && (T % DT) == 0) {
							if (chemoDosesGiven < 7) {
								ChemoTime[chemoDosesGiven] = T;
								chemoDosesGiven++;
							}
						}
					}
				}
			}


			const double R0 = CHEMO_RADIUS_CELLS;
			const double dt_hr = 1.0 / (60.0 * (double)DT);

			const double kmax_per_hr = 1.0 / 12.0;
			for (j = 0; j < chemoDosesGiven; ++j) {
				double t_hr = steps_to_hours(T - ChemoTime[j], DT);

				double Thalf_hours_j = (double)ThalfPerDose[j];
				if (Thalf_hours_j <= 0.0) { ChemoC[j] = 0.0; continue; }

				ChemoC[j] = Cblood_TMQ(t_hr, Cmax, Thalf_hours_j, Tmax_hours);
				if (!isfinite(ChemoC[j]) || ChemoC[j] < 0.0) ChemoC[j] = 0.0;

			}

			double cCell = 0.0;
			if (chemoDosesGiven > 0) {

				for (i = 0; i < numCells; ++i) {

					cCell = 0.0;
					for (j = 0; j < chemoDosesGiven; ++j) {
						double t_hr = steps_to_hours(T - ChemoTime[j], DT);
						if (t_hr < 0.0) continue;

						// distance from vessel wall (cellDist already subtracts Rvessel)
						double d = cellDist[i];
						if (d < 0.0) d = 0.0;

						//linear-to-zero at R0
						double w = 1.0 - (d / R0);
						if (w < 0.0) w = 0.0;

						cCell += ChemoC[j] * w;
					}
					cellConc[i] = cCell;

					double nHill = 1.0;
					double x = pow(cellConc[i], nHill);
					double y = pow(IC50,       nHill);
					double effect = (x + y > 0.0) ? (x / (x + y)) : 0.0;

					// Pharmacodynamic susceptibility window for TMZ-induced DNA lesion processing.
					double pd_window = 0.0;
					for (int jj = 0; jj < chemoDosesGiven; ++jj) {
						double dt_chemo_hr = steps_to_hours(T - ChemoTime[jj], DT);
						if (dt_chemo_hr >= 0.0) {
							double w = exp(-dt_chemo_hr / 24.0);
							if (w > pd_window) pd_window = w;
						}
					}

					effect *= pd_window;

					double cycling_susc =  (cellQ[i] == 1) ? 0.1 : 1.0;
					effect *= cycling_susc;

					double hazard_per_hr = kmax_per_hr * effect;
					if (hazard_per_hr < 0.0) hazard_per_hr = 0.0;

					double pkill = 1.0 - exp(-hazard_per_hr * dt_hr);
					if (pkill < 0.0) pkill = 0.0;
					if (pkill > 1.0) pkill = 1.0;

					if (randZerotoOne() <= pkill) {
						cellsToKill[numKillCells++] = i;
						Ch_rm++;
						cellType[i] = 3;
					}

				}
			}


			if (numKillCells > 0) {
				remove_dead_cells(hasExitedQ, &numCells,
						cellsX, cellsY,
						cellDeathTime, cellDivTime,
						cellConc,
						cellDiffTime, celldeDiffTime,
						cellType, cellDist,
						cellAge,  cellRadii,
						cellQ, radQ, radQt,
						cellsToKill, numKillCells);
			}


			/*****Handle Radiotherapy******/
			numKillCells = 0;
			giveRad = 0; 	
			double t_hours_event = steps_to_hours(T, DT);
			if (week < RadWeeks) {
				for (i=0; i< radDoses; i++) {
					if (day == radDay[i]) {
						if (hour == radHour[i]) {
							if (minute_of_hour == 0 && (T % DT) == 0) {
								if (radDoseA[i] > 0) { 
									radDose = radDoseA[i];
									giveRad = 1;
								}
							}
						}
					}
				}
			}
			if (giveRad == 1) {
				expected_kills=0.0;
				stemN =0;
				tbN = 0;
				int prevRadT = lastRadT; 
				TOT_RAD = TOT_RAD + 1;
				lastRadT = T;
				numKillCells=0;
				for (i=0; i < numCells; i++) {
					double alpha_tau, beta_tau, rho_tau;

					if (cellType[i] == 0) {             // GSC / stem-like
						alpha_tau = alphaS;
						beta_tau  = betaS;
						rho_tau   = RadRho;
						stemN++;
					} else {                            // tumor bulk / DTC
						alpha_tau = alphaT;
						beta_tau  = betaT;
						rho_tau   = 1.0;
						tbN++;
					}
					double C = cellConc[i];  // TMZ concentration at this cell at RT time (must match ksens units)
					double sens_factor = 2.0;
					if (sens_on && C > 0.0) {
						double s = C / (IC50 + C);
						// Time-local TMZ radiosensitization window.
						double recent_weight = 0.0;
						for (int jj = 0; jj < chemoDosesGiven; ++jj) {
							double dt_chemo_hr = steps_to_hours(T - ChemoTime[jj], DT);
							if (dt_chemo_hr >= 0.0) {
								double w = exp(-dt_chemo_hr / 24.0);
								if (w > recent_weight) recent_weight = w;
							}
						}

						sens_factor += 0.3 * s * recent_weight;
					}
					double expo_linear =  (alpha_tau * rho_tau * radDose) * sens_factor;
					double expo_quad   = (beta_tau  * rho_tau * radDose * radDose);
					double expo_total  = expo_linear + expo_quad;
					double pkill = 1.0 - exp(-expo_total);

					expected_kills += pkill;


					//using linear quadratic model of radiotherapy, determine likelihood of cell death
					if (cellType[i] == 0) {
						int fixed = (int) llround((60.0 * (double)DT) * LStem);     // QT_min in steps
						radQ[i]  = fixed;
						radQt[i] = 1;
					} else {
						int fixed = (int) llround((60.0 * (double)DT) * LTumorBulk);
						radQ[i]  = fixed;
						radQt[i] = 1;
					}
					if (randZerotoOne() <= pkill) {
						cellsToKill[numKillCells] = i;
						numKillCells = numKillCells + 1;
						cellType[i] = 3;
						RAD_KILL=RAD_KILL+1;
					}
					else {
						if (timeDepGamma == 1) {
							int t0 = (prevRadT < 0) ? -1 : (T - prevRadT);   //int t0 = T - lastRadT;
							double gamma0 = 0.4; // baseline fraction capable of reversion (eta0 / gamma0)
							if (t0 < 0) {
								gamma = gamma0; // first dose (no previous dose): use baseline
							} else {
								double dt = (double)t0 - mu;          // mu is in timesteps
								gamma = gamma0 * exp(-(dt*dt)/(2.0*sigmaS));
							}
						}

						celldeDiffTime[i] = -1;
						if (randZerotoOne() <= gamma) {
							if (cellQ[i] == 0) {
								if (cellType[i] < Zrevert) {
									if (cellType[i] == 1 ) celldeDiffTime[i] = exp_wait_steps(1.0 / aDeDiff, DT); 
								}
							}
						}
					}
				}
			}

			if (numKillCells > 0) {
				remove_dead_cells(hasExitedQ, &numCells, cellsX, cellsY, cellDeathTime, cellDivTime, cellConc, cellDiffTime, celldeDiffTime, cellType, cellDist, cellAge, cellRadii, cellQ, radQ, radQt, cellsToKill, numKillCells);
			}

			for (i=0; i < numCells; i++) {

				if (cellType[i] >1) cellType[i] = 1;
			}
			/*****Handle Quiescence******/
			for (i=0; i < numCells; i++) {
				if (radQ[i] > 0) {
					radQ[i] = radQ[i] - 1;
					radQt[i] = radQt[i] + 1;
					cellQ[i] = 1;
					//if leaving quiescence on next time step, set cell division rate to proliferate based on rTS and rTB
					if (radQ[i] <=0) {   
						hasExitedQ[i] = 1;

						double rate_i  = (cellType[i] == 0) ? rS : rTB;                 // /hour
						double mean_i  = (hasExitedQ[i] ? ((rate_i > 0.0) ? (1.0/rate_i) : 0.0)
								: 36.0);
						cellDivTime[i] = exp_wait_steps(mean_i, DT);
					}
				} else {
					radQt[i] = 0;
					cellQ[i] = 0;
				}
			}





			/*****Handle Stem Cell Retention Near Vessel Boundary******/
			for (i=0; i < numCells; i++) {
				if (cellDist[i] < (tau * R)) {
					if (cellType[i] == 1 && randZerotoOne() < 0.10) {
						cellType[i] = 0;
					}
				}			
			}
			/*****Handle Aging******/
			for (i=0; i < numCells; i++) {
				cellDeathTime[i] = cellDeathTime[i] - 1;
				if (cellQ[i] == 0)  {
					cellDivTime[i] = cellDivTime[i] - 1;
					if (cellDiffTime[i] > 0) cellDiffTime[i] = cellDiffTime[i] - 1;
					if (celldeDiffTime[i] > 0) celldeDiffTime[i] = celldeDiffTime[i] - 1;
				}

				if (cellDivTime[i] < 0) cellDivTime[i] = 0;
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

			age1Counter=0;
			dist1Counter=0;
			age2Counter=0;
			dist2Counter=0;
			typeCounter=0;
			maxDist = 0;
			for (i=0; i < numCells; i++)  {
				if (cellType[i] == 0) {
					age1Counter += cellAge[i]-cellDeathTime[i];
					dist1Counter += sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])-Rvessel; 
					if ( cellQ[i] >0) age1Counter += radQt[i];
				}
				else {
					age2Counter += cellAge[i]-cellDeathTime[i];
					dist2Counter += sqrt(cellsX[i]*cellsX[i] + cellsY[i] * cellsY[i])-Rvessel; 
					if ( cellQ[i] >0) age2Counter += radQt[i];
				}
				typeCounter += cellType[i];
				if (fabs(maxDist) < fabs(cellDist[i])) maxDist = cellDist[i];
			}
			if (T % 50 == 0) {
				fprintf(fp3, "%d, %d\n", T, numCells);
				fprintf(fp1, "%d,%d,%d,%d,%f,%d,%f,%f\n",
						T,
						numCells,
						stemCount[T],
						age1Counter,
						dist1Counter,
						age2Counter,
						dist2Counter,
						maxDist
				       );
			}	

		}//timesteps
		fclose(fp3);
		fclose(fp1);
	}//cycle
	fclose(fp);

	free(cellsToKill);
}

