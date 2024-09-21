/*Function to read in collision-induce absorption cross section data */
/* output opacity in cm+5 */
/* Interpolation to specified temperature in each layer */

#include <math.h>
#include "constant.h"

void readcia();

void readcia()
{	
	int i, j, s, nl;
	FILE *fim, *fout;
	char *temp1;
	double **cia;
	char crossfile[1024];
    double tfitting[zbin+1];
    for (i=1; i<=zbin; i++) {
        tfitting[i]=tl[i];
        if (tfitting[i]>2000.0) {
            tfitting[i]=2000.0;
        }
        if (tfitting[i]<100.0) {
            tfitting[i]=100.0;
        }
    }
	
	strcpy(crossfile,CROSSHEADING);
	strcat(crossfile,"H2-H2_CIA.dat");
	printf("%s\n",crossfile);
	fim = fopen(crossfile,"r");
	s=LineNumber(fim, 10000);
	fclose(fim);
	nl = s-1;
    double wave[nl];
    double cross[nl];
    double temp[20] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0};
    char   dataline[10000];
    
    cia = dmatrix(0,nl-1,0,19);
	fim = fopen(crossfile,"r");
	temp1=fgets(dataline, 10000, fim); /* Read in the header line */
	for (i=0; i<nl; i++) {
		fscanf(fim, "%lf", wave+i);
		for (j=0; j<19; j++) {
			fscanf(fim, "%le", &cia[i][j]);
		}
		fscanf(fim, "%le\n", &cia[i][19]);
	}
	fclose(fim);
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
			}
	}
	fout=fopen("AuxillaryOut/CheckCIA_H2H2.dat","w");
	for (j=0; j<NLAMBDA; j++) {
		fprintf(fout,"%2.6f\t",wavelength[j]);
		for (i=1; i<=zbin; i++) {
			fprintf(fout,"%2.6e\t",H2H2CIA[i][j]);
		}
		fprintf(fout,"\n");
	}
	fclose(fout);
	free_dmatrix(cia,0,nl-1,0,19);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"H2-He_CIA.dat");
    printf("%s\n",crossfile);
    cia = dmatrix(0,nl-1,0,19);
    fim = fopen(crossfile,"r");
    temp1=fgets(dataline, 10000, fim); /* Read in the header line */
    for (i=0; i<nl; i++) {
        fscanf(fim, "%lf", wave+i);
        for (j=0; j<19; j++) {
            fscanf(fim, "%le", &cia[i][j]);
        }
        fscanf(fim, "%le\n", &cia[i][19]);
    }
    fclose(fim);
    /* Calculate the cross sections */
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j++) {
            H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
        }
    }
    fout=fopen("AuxillaryOut/CheckCIA_H2He.dat","w");
    for (j=0; j<NLAMBDA; j++) {
        fprintf(fout,"%2.6f\t",wavelength[j]);
        for (i=1; i<=zbin; i++) {
            fprintf(fout,"%2.6e\t",H2HeCIA[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    free_dmatrix(cia,0,nl-1,0,19);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"H2-H_CIA.dat");
    printf("%s\n",crossfile);
    cia = dmatrix(0,nl-1,0,19);
    fim = fopen(crossfile,"r");
    temp1=fgets(dataline, 10000, fim); /* Read in the header line */
    for (i=0; i<nl; i++) {
        fscanf(fim, "%lf", wave+i);
        for (j=0; j<19; j++) {
            fscanf(fim, "%le", &cia[i][j]);
        }
        fscanf(fim, "%le\n", &cia[i][19]);
    }
    fclose(fim);
    /* Calculate the cross sections */
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j++) {
            H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
        }
    }
    fout=fopen("AuxillaryOut/CheckCIA_H2H.dat","w");
    for (j=0; j<NLAMBDA; j++) {
        fprintf(fout,"%2.6f\t",wavelength[j]);
        for (i=1; i<=zbin; i++) {
            fprintf(fout,"%2.6e\t",H2HCIA[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    free_dmatrix(cia,0,nl-1,0,19);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"N2-N2_CIA.dat");
    printf("%s\n",crossfile);
    cia = dmatrix(0,nl-1,0,19);
    fim = fopen(crossfile,"r");
    temp1=fgets(dataline, 10000, fim); /* Read in the header line */
    for (i=0; i<nl; i++) {
        fscanf(fim, "%lf", wave+i);
        for (j=0; j<19; j++) {
            fscanf(fim, "%le", &cia[i][j]);
        }
        fscanf(fim, "%le\n", &cia[i][19]);
    }
    fclose(fim);
    /* Calculate the cross sections */
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j++) {
            N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
        }
    }
    fout=fopen("AuxillaryOut/CheckCIA_N2N2.dat","w");
    for (j=0; j<NLAMBDA; j++) {
        fprintf(fout,"%2.6f\t",wavelength[j]);
        for (i=1; i<=zbin; i++) {
            fprintf(fout,"%2.6e\t",N2N2CIA[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    free_dmatrix(cia,0,nl-1,0,19);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"N2-H2_CIA.dat");
    printf("%s\n",crossfile);
    cia = dmatrix(0,nl-1,0,19);
    fim = fopen(crossfile,"r");
    temp1=fgets(dataline, 10000, fim); /* Read in the header line */
    for (i=0; i<nl; i++) {
        fscanf(fim, "%lf", wave+i);
        for (j=0; j<19; j++) {
            fscanf(fim, "%le", &cia[i][j]);
        }
        fscanf(fim, "%le\n", &cia[i][19]);
    }
    fclose(fim);
    /* Calculate the cross sections */
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j++) {
            N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
        }
    }
    fout=fopen("AuxillaryOut/CheckCIA_N2H2.dat","w");
    for (j=0; j<NLAMBDA; j++) {
        fprintf(fout,"%2.6f\t",wavelength[j]);
        for (i=1; i<=zbin; i++) {
            fprintf(fout,"%2.6e\t",N2H2CIA[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    free_dmatrix(cia,0,nl-1,0,19);
    
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"CO2-CO2_CIA.dat");
    printf("%s\n",crossfile);
    cia = dmatrix(0,nl-1,0,19);
    fim = fopen(crossfile,"r");
    temp1=fgets(dataline, 10000, fim); /* Read in the header line */
    for (i=0; i<nl; i++) {
        fscanf(fim, "%lf", wave+i);
        for (j=0; j<19; j++) {
            fscanf(fim, "%le", &cia[i][j]);
        }
        fscanf(fim, "%le\n", &cia[i][19]);
    }
    fclose(fim);
    /* Calculate the cross sections */
    for (i=1; i<=zbin; i++) {
        for (j=0; j<NLAMBDA; j++) {
            CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], wave, nl, temp, 20, cia);
        }
    }
    fout=fopen("AuxillaryOut/CheckCIA_CO2CO2.dat","w");
    for (j=0; j<NLAMBDA; j++) {
        fprintf(fout,"%2.6f\t",wavelength[j]);
        for (i=1; i<=zbin; i++) {
            fprintf(fout,"%2.6e\t",CO2CO2CIA[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    free_dmatrix(cia,0,nl-1,0,19);
	
}
