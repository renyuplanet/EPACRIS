/*Function to read in cross section data */
/* input cross section in m^2 */
/* input wavelength in microns */
/* output mean cross section in cm^2 */
/* Interpolation to specified T and P */

#include <math.h>
#include "constant.h"

void readcross(char Fname[], double **xsc);

void readcross(char Fname[], double **xsc)
{
	int i, j, k;
	FILE *fim;
	double ***opac;
	opac = f3tensor(0,NLAMBDA-1,0,NTEMP-1,0,NPRESSURE-1);
	double wave[NLAMBDA]; /* in m */
	double wave1[NLAMBDA]; /* in nm */
	double temp[NTEMP]; /* in k */
	double pres[NPRESSURE]; /* in Pa */
	double presdummy, cross[NLAMBDA];
	double pfitting, tfitting;
	
//ms2022	printf("%s\n", Fname);
	fim = fopen(Fname,"r");	
	/* Header Lines */
	for (i=0; i<NTEMP-1; i++) {
		fscanf(fim, "%lf", temp+i);
	}
	i=NTEMP-1;
	fscanf(fim, "%lf\n", temp+i);
	for (i=0; i<NPRESSURE-1; i++) {
		fscanf(fim, "%le", pres+i);
	}
	i=NPRESSURE-1;
	fscanf(fim, "%le\n", pres+i);
	/* Read in data */
	for (i=0; i<NLAMBDA; i++) {
		fscanf(fim, "%le\n", wave+i);
		for (j=0; j<NPRESSURE; j++) {
			fscanf(fim, "%le", &presdummy);
			for (k=0; k<NTEMP-1; k++) {
				fscanf(fim, "%le", &opac[i][k][j]);
			}
			k=NTEMP-1;
			fscanf(fim, "%le\n", &opac[i][k][j]);
		}
	}
	fclose(fim);
	
	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = wave[i]*1.0E+9; /* convert to nm */
    }
    
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i]>pres[NPRESSURE-1]) {
				pfitting=pres[NPRESSURE-1];
            } else if (pl[i]<pres[0]) {
                pfitting=pres[0];
            } else {
                pfitting=pl[i];
            }
			if (tl[i]>temp[NTEMP-1]) {
				tfitting=temp[NTEMP-1];
			} else if (tl[i]<temp[0]) {
                tfitting=temp[0];
            } else {
				tfitting=tl[i];
			}
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, temp, NTEMP, pres, NPRESSURE, *(opac+j))*1.0E+4,0.0); /* convert to cm^2 */
            xsc[i][j] = cross[j];
		}
		/*Interpolation(wavelength, NLAMBDA, *(xsc+i), wave1, cross, NLAMBDA, 0);*/
	}
    /*for (i=0; i<NLAMBDA; i++) { printf("%e %e\n", wavelength[i], xsc[zbin][i]);}*/

	/* free_dmatrix(opacj,1,NTEMP,1,NPRESSURE);*/ 
	free_f3tensor(opac,0,NLAMBDA-1,0,NTEMP-1,0,NPRESSURE-1);

}
