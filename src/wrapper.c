/*Copyright (c) 2018, Dana Farber Cancer Institute, Duke University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the DFCI or Duke nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL DANA FARBER CANCER INSTITUTE OR DUKE UNIVERSITY
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NO T LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "math.h"
#include "time.h"
#include "cellModelLoop.h"


//Note: Run this with passing in the fileName for the schedule
int main(int argc, char** argv) {


	/****************** Open Input File ******************************/
	FILE *fp,*fpsched,*fout;
	fp =fopen(argv[1],"r");  //open file that was passed in as argv[1]
	if (fp==NULL){      //error checking
		printf("Error! opening input file");
		exit(1);
	}

	/****************** Open Output File ******************************/
		fout =fopen("full_output.csv","w");  //open file that was passed in as argv[1]
		if (fp==NULL){      //error checking
			printf("Error! opening input file");
			exit(1);
		}

	/*************** Read In Input File Parameters *******************/
	char name[256];
	char * token;
	char fileName[256];
	char outputFileName[256];
	int maxTimeSteps;
	int Z, T;
	int Zrevert;

	fscanf(fp,"%s = %s\n",name,fileName);
	fscanf(fp,"%s = %d\n",name,&maxTimeSteps);
	fscanf(fp,"%s = %d\n",name,&Z);
	fscanf(fp,"%s = %d\n",name,&Zrevert);
	fclose(fp); //close input file

	//print input and schedule files being used to stdout
	printf("Input file: %s\n", argv[1]);
	printf("Schedule file: %s\n", fileName);

	/*************** Read In Schedule Parameters *******************/
	fp =fopen(fileName,"r");  //open file that was passed in as the schedule
	if (fp==NULL){      //error checking
		printf("Error! opening schedule file");
		exit(1);
	}
	int ChemoWeeks;
	int ChemoDays;
	int ChemoDoses = ChemoDays;   //note: hard-coded!
	int ThalfPerDose[7];
	int ChemoDay[7];
	int ChemoHour[7];
	double ChemoC[7];
	double Cmax;
	int RadWeeks;
	int radDoses;
	int radDose;
	int radDoseA[10];
	int radHour[10];
	int radDay[10];

	fscanf(fp,"%s = %d\n",name,&ChemoWeeks);  //start of chemo parameters
	fscanf(fp,"%s = %d\n",name,&ChemoDays);
	fscanf(fp,"%s = %d %d %d %d %d %d %d\n",name,&ThalfPerDose[0],&ThalfPerDose[1],&ThalfPerDose[2],&ThalfPerDose[3],&ThalfPerDose[4],&ThalfPerDose[5],&ThalfPerDose[6]);
	fscanf(fp,"%s = %d %d %d %d %d %d %d\n",name,&ChemoDay[0],&ChemoDay[1],&ChemoDay[2],&ChemoDay[3],&ChemoDay[4],&ChemoDay[5],&ChemoDay[6]);
	fscanf(fp,"%s = %d %d %d %d %d %d %d\n",name,&ChemoHour[0],&ChemoHour[1],&ChemoHour[2],&ChemoHour[3],&ChemoHour[4],&ChemoHour[5],&ChemoHour[6]);
	fscanf(fp,"%s = %lf %lf %lf %lf %lf %lf %lf\n",name,&ChemoC[0],&ChemoC[1],&ChemoC[2],&ChemoC[3],&ChemoC[4],&ChemoC[5],&ChemoC[6]);
	fscanf(fp,"%s = %lf\n",name,&Cmax);
	fscanf(fp,"%s = %d\n",name,&RadWeeks);  //start of radation parameters
	fscanf(fp,"%s = %d\n",name,&radDoses);
	fscanf(fp,"%s = %d\n",name,&radDose);
	fscanf(fp,"%s = %d %d %d %d %d %d %d %d %d %d\n",name,&radDoseA[0],&radDoseA[1],&radDoseA[2],&radDoseA[3],&radDoseA[4],&radDoseA[5],&radDoseA[6],&radDoseA[7],&radDoseA[8],&radDoseA[9]);
	fscanf(fp,"%s = %d %d %d %d %d %d %d %d %d %d\n",name,&radHour[0],&radHour[1],&radHour[2],&radHour[3],&radHour[4],&radHour[5],&radHour[6],&radHour[7],&radHour[8],&radHour[9]);
	fscanf(fp,"%s = %d %d %d %d %d %d %d %d %d %d\n",name,&radDay[0],&radDay[1],&radDay[2],&radDay[3],&radDay[4],&radDay[5],&radDay[6],&radDay[7],&radDay[8],&radDay[9]);
	int i;

	outputFileName[0] = 0;
	strcat(outputFileName,"output");
	strcat(outputFileName + strlen(outputFileName),argv[1]+5);

	int radDoseAg[10];
	int radHourg[10];
	int radDayg[10];
	int radDosesg=radDoses;

	for (i=0; i<10; i++) {
		radDoseAg[i] = radDoseA[i];
		radHourg[i] = radHour[i];
		radDayg[i] = radDay[i];
	}

	int *tbCount; 
	tbCount = (int *) malloc(sizeof(int)*maxTimeSteps);
	int *stemCount; //array to keep track of the number of stem cells over time
	stemCount = (int *) malloc(sizeof(int)*maxTimeSteps);

	//Initialize counts
	for (T=0; T < maxTimeSteps; T++) {
		tbCount[T] = 0;
		stemCount[T] = 0;
	}


	cellModelLoop(maxTimeSteps,Z, Zrevert,ChemoWeeks,ChemoDays, ChemoDoses, ThalfPerDose, ChemoDay, ChemoHour, ChemoC, Cmax, RadWeeks, radDoses, radDose, radDoseA, radHour, radDay, stemCount,tbCount);


	for (i=0; i < maxTimeSteps; i++) {
		fprintf(fout,"%d, %d, %d\n", i, (int) floor(tbCount[i]),  (int) floor(stemCount[i]) );
	}
	fclose(fout); //close input file
	free(tbCount);
	free(stemCount);
}
