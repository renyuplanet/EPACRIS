//
// Original Author: 
char authorig[] = "Renyu Hu (renyu.hu@jpl.nasa.gov)";
// Version: 
char version[] = "2023, Major RadTrans and Convection upgrades";
// Editor v2023: 
char editor[] = "Markus Scheucher (markus.scheucher@jpl.nasa.gov)";
/* Modifier symbols: 
    * 'ms_' ... in file names
    * 'markus<date>' ... in description lines/blocks 
    * 'ms<date>' ... in-line comments
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//#include "Input/debugging.h"
//
//#include "Input/jupiter-guillot.h"
//#include "Input/jupiter-time.h"
//#include "Input/jupiter-jacob.h"
//
//#include "Input/k218b-guillot.h"
//#include "Input/k218b-time.h"
//#include "Input/k218b-jacob.h"
//
#include "Input/toi1634b-jacob-200.h"

#include "constant.h"
#include "ms_functions.h" //ms2021
#include "routine.h"
#include "global_temp.h"
#include "nrutil.h"

double TAUdoub[2*zbin+1], Tdoub[2*zbin+1], Pdoub[2*zbin+1], MMdoub[2*zbin+1], zdoub[2*zbin+1]; //double grid for non-isothermal layers: all bottom-up!!
double meanmolecular[zbin+1];
double zl[zbin+1];
double pl[zbin+1];
double tl[zbin+1];
double MM[zbin+1];
double MMZ[zbin+1];
double wavelength[NLAMBDA];
double solar[NLAMBDA];
double crossr[NLAMBDA], crossa[3][NLAMBDA], sinab[3][NLAMBDA], asym[3][NLAMBDA];
double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO;
double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4;
double **opacNH3;
double **opacC2H2, **opacC2H4, **opacC2H6, **opacHCN, **opacCH2O2, **opacHNO3;
double **opacN2O, **opacN2, **opacNO, **opacNO2, **opacOCS;
double **opacHF, **opacHCl, **opacHBr, **opacHI, **opacClO, **opacHClO;
double **opacHBrO, **opacPH3, **opacCH3Cl, **opacCH3Br, **opacDMS, **opacCS2;
double **cross, **crosst; //PhotoCross in Clima
int *stdcross; //PhotoCross in Clima
double *qysum; //PhotoCross in Clima
double MeanCO2[zbin+1], MeanO2[zbin+1], MeanSO2[zbin+1], MeanH2O[zbin+1], MeanOH[zbin+1], MeanH2CO[zbin+1];
double MeanH2O2[zbin+1], MeanHO2[zbin+1], MeanH2S[zbin+1], MeanCO[zbin+1], MeanO3[zbin+1], MeanCH4[zbin+1];
double MeanNH3[zbin+1];	
double MeanC2H2[zbin+1], MeanC2H4[zbin+1], MeanC2H6[zbin+1], MeanHCN[zbin+1], MeanCH2O2[zbin+1], MeanHNO3[zbin+1];
double MeanN2O[zbin+1], MeanN2[zbin+1], MeanNO[zbin+1], MeanNO2[zbin+1], MeanOCS[zbin+1];

double SMeanCO2[zbin+1], SMeanO2[zbin+1], SMeanSO2[zbin+1], SMeanH2O[zbin+1], SMeanOH[zbin+1], SMeanH2CO[zbin+1];
double SMeanH2O2[zbin+1], SMeanHO2[zbin+1], SMeanH2S[zbin+1], SMeanCO[zbin+1], SMeanO3[zbin+1], SMeanCH4[zbin+1];
double SMeanNH3[zbin+1];
double SMeanC2H2[zbin+1], SMeanC2H4[zbin+1], SMeanC2H6[zbin+1], SMeanCH2O2[zbin+1];
double SMeanHCN[zbin+1], SMeanN2O[zbin+1], SMeanN2[zbin+1], SMeanNO[zbin+1], SMeanNO2[zbin+1], SMeanOCS[zbin+1], SMeanHNO3[zbin+1];

int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
int    numr=0, numm=0, numt=0, nump=0, numx=0, numc=0, numf=0, numa=0, waternum=0, waterx=0;
double xx[zbin+1][NSP+1];
double clouds[zbin+1][NSP+1] = {0.0}; //ms22: moist adiabat & clouds update
double mkv[zbin+1], Tnew[zbin+1], Pnew[zbin+1];

double H2H2CIA[zbin+1][NLAMBDA], H2HeCIA[zbin+1][NLAMBDA], H2HCIA[zbin+1][NLAMBDA], N2H2CIA[zbin+1][NLAMBDA], N2N2CIA[zbin+1][NLAMBDA], CO2CO2CIA[zbin+1][NLAMBDA];
double MeanH2H2CIA[zbin+1], MeanH2HeCIA[zbin+1], MeanH2HCIA[zbin+1], MeanN2H2CIA[zbin+1], MeanN2N2CIA[zbin+1],MeanCO2CO2CIA[zbin+1];
double SMeanH2H2CIA[zbin+1], SMeanH2HeCIA[zbin+1], SMeanH2HCIA[zbin+1], SMeanN2H2CIA[zbin+1], SMeanN2N2CIA[zbin+1], SMeanCO2CO2CIA[zbin+1];
/*markus2021
char filleq[] = "======================================================================";
char fillmi[] = "----------------------------------------------------------------------";
char fillst[] = "**********************************************************************";
char fillpl[] = "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
*/
//end edits
double rt_drfluxmax_init;
int RTstepcount;
double THETAREF; //previously in input file, but redefined there as 'degrees' rather than radian

#include "GetData.c"
#include "chemequil.c"
#include "Interpolation.c"
#include "nrutil.c"
#include "Convert.c"
#include "TPPara.c"
//#include "gasheat.c"
#include "ms_gasheat.c"
#include "GreyTemp.c"
#include "ms_Climate_test_v1.c"

#include "RefIdx.c"
#include "planckmean.c"
#include "planckmeanCIA.c"
#include "readcia.c"
#include "readcross.c"
#include "printout_std.c"
#include "printout_std_t_exp.c"

//=== START MAIN PROGRAM =================================
//========================================================
void main(int argc, char *argv[]) //ms2022: getting rid of warnings
{
//========================================================
//========================================================
    double time,t_passed,t2_passed;
    struct timespec tstart,tend;
    //clock_gettime(CLOCK_REALTIME, &tstart);
    //clock_gettime(CLOCK_REALTIME, &tend);
    //t_passed = ((double)tend.tv_sec*1e9 + tend.tv_nsec) - ((double)tstart.tv_sec*1e9 + tstart.tv_nsec);
    //printf("\n%s %f\n","task took [ms]:", t_passed/1e6); 

    THETAREF = PI/180.0*THETAANGLE; // comparing with original
    //THETAREF = 1.0471;
//markus2021 
    printf("%s\n",fillst);
    printf("%s", "ExoSurfAir launched on ");
    datetime();
    printf("%s %s\n", "Original Author: ", authorig);
    printf("%s %s\n", "Version: ", version);
    printf("%s %s\n", "Edited by: ", editor);
    printf("%s\n\n",fillst);
   
    if (strlen(COMMENTS)>0){
    printf("%s %s\n", "Special Comments: ",COMMENTS);
    printf("%s\n\n",fillmi);
    }
    
    printf("%s %s\n", "Input File: ",IN_FILE_NAME);
    printf("%s %s\n", "Results Directory: ",OUT_DIR);
    printf("%s\n\n",fillmi);

    printf("%s\n", "PLANETARY SETUP:");
    printf("%s %s\n", "Host Star = ", STAR_SPEC);
    printf("%s %e %s\n", "M_planet = ", MASS_PLANET, "kg" );
    printf("%s %e %s\n", "Interior Temperature = ", TINTSET, "K" );
    printf("%s %f\n", "Advection Factor = ", FADV);
    printf("%s %f %s%f%s\n", "Zenith Angle = ", THETAANGLE, "(",THETAREF," rad)");

    printf("%s\n\n",fillmi);
    
    printf("%s\n", "CLIMATE SETUP:");
    if (TIME_STEPPING == 0) printf("%s %f %s %e %s\n", "Equilibrium solver using flux Jacobian with delta_Temp=", RJACOB_TEMPVAR,"& max step = ",R_RELAX,"Temp_old");
    if (TIME_STEPPING == 1) printf("%s\n", "Time stepping towards equilibrium");
    if (TWO_STR_SOLVER == 1) 
    {
        printf("%s\n","Using Heng+2018 2-STREAM equations for anisotropic scattering and non-isothermal layers (by M.Scheucher)");
        if (RT_FLUX_SOLVER <= 2) printf("%s\n","Equations 25&26 re-written into tri-diagonal form (see ms_2stream.c)");
        if (RT_FLUX_SOLVER >= 3) printf("%s\n","Equations 25&26 used in original penta-diagonal form (see ms_2stream.c)");
        printf("%s%d%s","RT_FLUX_SOLVER of choice: [",RT_FLUX_SOLVER,"]");
        if (RT_FLUX_SOLVER == 0) printf("%s"," : Thomas algorithm (BOA->TOA & back-substitution TOA->BOA)");
        if (RT_FLUX_SOLVER == 1) printf("%s"," : Thomas algorithm (TOA->BOA & back-substitution BOA->TOA)");
        if (RT_FLUX_SOLVER == 2) printf("%s"," : LU decomposition (ludcmp&lubksb)");
        if (RT_FLUX_SOLVER == 3) printf("%s"," : Block Tridiagonal solver using LU decomposition on 2x2 blocks");
        if (RT_FLUX_SOLVER == 4) printf("%s"," : PTRANS-I pentadiagonal solver (Askar & Karawia, 2015)");
        if (RT_FLUX_SOLVER == 5) printf("%s"," : Sogabe 2008 pentadiagonal solver (Algorithms 2&3)");
        printf("\n");
    }
    if (TWO_STR_SOLVER == 0) printf("%s\n","Using Toon+89 Delta-2-STREAM Solver (by R.Hu)");
    if (CONVEC_ADJUST == 1) printf("%s\n","Convective Adjustment: Conserving Enthalpy (using potential temperature) (by M.Scheucher)");
    if (CONVEC_ADJUST == 0) printf("%s\n","Convective Adjustment: Dry with simple R/cp (by R.Hu)");
    printf("%s\n\n",fillmi);
//end edits
//atexit(pexit);exit(0); //ms debugging mode
	  int s,i,ii,j,jj,jjj,k,nn,qytype,stdnum,iradmax;
	  int numr1=1,numm1=1,nump1=1,numt1=1;
	  int nums;
	  int numx1=1, numf1=1, numc1=1;
	  char *temp;
	  char dataline[10000];
	  double temp1, wavetemp, crosstemp, DD, GA, DenZ;
	  double z[zbin+1], T[zbin+1], PP[zbin+1], P[zbin+1], tempeq[zbin+1];
//          double **JJ, **cross, **qy, *wavep, *crossp, *qyp, *qyp1, *qyp2, *qyp3, *qyp4, *qyp5, *qyp6, *qyp7;
//	  double **crosst, **qyt, *crosspt, *qypt, *qyp1t, *qyp2t, *qyp3t, *qyp4t, *qyp5t, *qyp6t, *qyp7t;
          double **JJ, **qy, *wavep, *crossp, *qyp, *qyp1, *qyp2, *qyp3, *qyp4, *qyp5, *qyp6, *qyp7;
	  double **qyt, *crosspt, *qypt, *qyp1t, *qyp2t, *qyp3t, *qyp4t, *qyp5t, *qyp6t, *qyp7t;
	  FILE *fspecies, *fzone, *fhenry, *fp, *fp1, *fp2, *fp3;
          FILE *fout, *fout1, *fout21, *fout22, *fout3, *fout4, *fcheck, *ftemp, *fout5, *foutp, *fcheckgibbs;
	  double **mixequil;

      GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET; /* Planet Surface Gravity Acceleration, in SI */	
	
	/*Set the wavelength for calculation*/
	double dlambda, start, interval, lam[NLAMBDA];
	start = log10(LAMBDALOW);
	interval = log10(LAMBDAHIGH) - log10(LAMBDALOW);
	dlambda = interval / (NLAMBDA-1.0);
	for (i=0; i<NLAMBDA; i++){
		wavelength[i] = pow(10.0, start+i*dlambda)*1.0E+9; /* in nm */
		lam[i] = wavelength[i]*1.0E-3; /* in microns */
                /*markus2021*/
                //if (i==4530 || i==4531) printf("%s\t%i\t%s%f\n","wavelength #",i,"= ",wavelength[i]);
	}
    printf("%s\n", "Set wavelength done");
	
	/* Rayleigh Scattering */
//ms2023: to be revisited regarding pressure dependencies
	double refidx0,DenS;
	DenS=101325.0/KBOLTZMANN/273.0*1.0E-6; 
	for (i=0; i<NLAMBDA; i++){
		if (RefIdxType == 0) { refidx0=AirRefIdx(wavelength[i]);}
		if (RefIdxType == 1) { refidx0=CO2RefIdx(wavelength[i]);}
		if (RefIdxType == 2) { refidx0=HeRefIdx(wavelength[i]);}
		if (RefIdxType == 3) { refidx0=N2RefIdx(wavelength[i]);}
		if (RefIdxType == 4) { refidx0=NH3RefIdx(wavelength[i]);}
		if (RefIdxType == 5) { refidx0=CH4RefIdx(wavelength[i]);}
		if (RefIdxType == 6) { refidx0=H2RefIdx(wavelength[i]);}
		if (RefIdxType == 7) { refidx0=O2RefIdx(wavelength[i]);}
		if (refidx0 < 1.0) { refidx0 = 1.0; }
		crossr[i]=1.061*8.0*pow(PI,3)*pow(pow(refidx0,2)-1,2)/3.0/pow(wavelength[i]*1.0E-7,4)/DenS/DenS;
        if (RefIdxType == 6) {crossr[i] = 8.14e-13*pow(wavelength[i]*10.0,-4)+1.28e-6*pow(wavelength[i]*10.0,-6)+1.61*pow(wavelength[i]*10.0,-8); } /* Dalgarno 1962 */
		/* printf("%s\t%f\t%s\t%e\n", "The reyleigh scattering cross-section at wavelength", wavelength[i], "nm is", crossr[i]); */
	}
    printf("%s\n", "Rayleigh scattering done");
	/* Obtain the stellar radiation */
	fp2=fopen(STAR_SPEC,"r");
	fp3=fopen(STAR_SPEC,"r");
	s=LineNumber(fp2, 1000);
	double swave[s], sflux[s];
	GetData(fp3, 1000, s, swave, sflux);
	fclose(fp2);
	fclose(fp3);
	Interpolation(wavelength, NLAMBDA, solar, swave, sflux, s, 0);
	for (i=0; i<NLAMBDA; i++) {
		solar[i] = solar[i]/ORBIT/ORBIT*FaintSun;  /* convert from flux at 1 AU */
	}
	i=0;
	while (solar[i]>0 || wavelength[i]<9990 ) { i++;}
	for (j=i; j<NLAMBDA; j++) {
		solar[j] = solar[i-1]*pow(wavelength[i-1],4)/pow(wavelength[j],4); //ms2021: what?
	}
	printf("%s\n", "Stellar radiation data is imported.");
	
	/* Initial Mean Molecular Mass */
	for (j=1; j<=zbin; j++) {
		meanmolecular[j] = AIRM;
	}
	
printf("%s\n\n",fillmi);
	double PMIN, PMAX, PSTEP, scaleheight;
	/* Set up the P-T-z for calculation */
	if (TPMODE==1) {
		fp=fopen(TPLIST,"r");
		fp1=fopen(TPLIST,"r");
		s=LineNumber(fp, 1000);
		double Height[s];
		double Temp[s];
		double Pre[s], lPre[s];
		GetData3(fp1, 1000, s, Height, Pre, Temp);
		fclose(fp);
		fclose(fp1);
		Reverse(Pre,s);
		Reverse(Temp,s);
		PMIN = Pre[0]; //input as P = pow(10,Pre)
		PMAX = Pre[s-1]; //input as P = pow(10,Pre)
//ms2023		PSTEP = (PMAX - PMIN)/zbin;
        for (j=0;j<s;j++) lPre[j] = log(pow(10,Pre[j]));
                PSTEP = log( pow(10,PMAX) / pow(10,PMIN) ) / zbin; //convert to natural log for obvious reasons
		for (j=0; j<=zbin; j++) { 
//ms2023			PP[j] = PMAX - j*PSTEP;
			PP[j] = log(pow(10,PMAX)) - j*PSTEP;
//ms2023			P[j]  = pow(10.0, PP[j]);
                        P[j] = exp(PP[j]); //natural log again
                        Pdoub[2*j] = P[j]; //ms2023: double grid
                        Tdoub[2*j] = T[j]; //ms2023: double grid
		}
//ms2023		Interpolation(PP,zbin+1,T,Pre,Temp,s,0);
		Interpolation(PP,zbin+1,T,lPre,Temp,s,2); //here interpolated linear in log(p) rather than 10^PP // rh2024: allow extrapolation to nearest
// rh2024        for (j=0; j<=zbin; j++) {printf("%s %f %f\n","PT",PP[j],T[j]);}
		for (j=1; j<=zbin; j++) {
			tl[j] = (T[j]+T[j-1])/2.0; /* Temperature at the center of layer */
//ms2023			pl[j] = sqrt(P[j]*P[j-1]); /* Pressure at the center of layer */
			pl[j] = exp( (PP[j]+PP[j-1]) /2.0 ); // center pressure is like double resolution grid, ergo half log(P)
                        
                        Pdoub[2*j-1] = pl[j]; //ms2023: double grid
                        Tdoub[2*j-1] = tl[j]; //ms2023: double grid
		}
		z[0] = 0.0;
                //zdoub[0] = z[0];
		for (j=1; j<=zbin; j++) {
			scaleheight = KBOLTZMANN * tl[j] / meanmolecular[j] / AMU / GA /1000.0 ; /* km */
			z[j] = z[j-1] - scaleheight*log(P[j]/P[j-1]);
			zl[j] = z[j-1] - scaleheight*log(pl[j]/P[j-1]); //ms2023
                        zdoub[2*j] = z[j]; //ms2023: double grid
		}
//ms2023		for (j=1; j<=zbin; j++) {
//ms2023			zl[j] = (z[j]+z[j-1])/2.0;
                        //zdoub[2*j-1] = zl[j]; //ms2023: double grid
//ms2023		}
                printf("%s %s\n", "Initial TP profile imported from ",TPLIST);
	}
	
	/* Set up the P-T-z for calculation */
	if (TPMODE==0) {
		TPPara(P,T,TINV,zbin+1,PTOP,TTOP,PMIDDLE,TMIDDLE,PSTR,TSTR,PTROP,TTROP,PBOTTOM);
		Reverse(P,zbin+1); //here P is in Pa
		Reverse(T,zbin+1);
        /*for (j=0; j<=zbin; j++) {
            printf("%le %lf\n", P[j], T[j]);
        }*/
		for (j=1; j<=zbin; j++) {
			tl[j] = (T[j]+T[j-1])/2.0; /* Temperature at the center of layer */
//ms2023			pl[j] = sqrt(P[j]*P[j-1]); /* Pressure at the center of layer */
			pl[j] = exp( log( P[j] * P[j-1] ) / 2.0 ); // center pressure is like double resolution grid, ergo half log(P)
                        Pdoub[2*j] = P[j]; //ms2023: double grid
                        Tdoub[2*j] = T[j]; //ms2023: double grid
                        Pdoub[2*j-1] = pl[j]; //ms2023: double grid
                        Tdoub[2*j-1] = tl[j]; //ms2023: double grid
		}
                Pdoub[0] = P[0]; //ms2023: double grid
                Tdoub[0] = T[0]; //ms2023: double grid
		z[0] = 0.0;
                //zdoub[0] = z[0];
		for (j=1; j<=zbin; j++) {
			scaleheight = KBOLTZMANN * tl[j] / AIRM / AMU / GA /1000.0 ; /* km */
			z[j] = z[j-1] - scaleheight*log(P[j]/P[j-1]);
			zl[j] = z[j-1] - scaleheight*log(pl[j]/P[j-1]); //ms2023
                        zdoub[2*j] = z[j]; //ms2023: double grid
		}
//ms2023		for (j=1; j<=zbin; j++) {
//ms2023			zl[j] = (z[j]+z[j-1])/2.0;
                        //zdoub[2*j-1] = zl[j]; //ms2023: double grid
//ms2023		}
                printf("%s\n%-10s  %-10s\n", "Initial TP profile calculated:","T[K]","P[Pa]");
		for (j=zbin; j>=0; j--) {printf("%-10.2f  %-10.2e\n",T[j],P[j]);}
        }

	FILE *TPPrint;
	TPPrint=fopen("AuxillaryOut/TPCheck.dat","w");
	for (j=0; j<=zbin; j++) {
		MMZ[j] = P[j]/KBOLTZMANN/T[j]*1.0E-6;
                MMdoub[2*j] = MMZ[j]; //ms2023: double grid
	}
        printf("%s\n%-10s %-10s %-10s %-14s\n", "Mid-layer TP profile:","Height[km]","P[Pa]","T[K]","Dens[mol/cm3]");
	for (j=zbin; j>=1; j--) {
		MM[j]=pl[j]/KBOLTZMANN/tl[j]*1.0E-6; /*unit: Molecule cm-3*/
                MMdoub[2*j-1] = MM[j]; //ms2023: double grid
        printf("%-10.1f %-10.2e %-10.2f %-14.2e\n", zl[j], pl[j], tl[j], MM[j]);
	fprintf(TPPrint, "%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
	} 
	printf("%s\n", "The Z-T-P data is imported/calculated.");
	fclose(TPPrint);
//atexit(pexit);exit(0); //ms debugging mode
	
printf("%s\n\n",fillmi);
	/* Get the species list */
	fspecies=fopen(SPECIES_LIST, "r");
	fout21=fopen(OUT_FILE1,"w");
	fout22=fopen(OUT_FILE2,"w");
    fprintf(fout21, "%s\t\t\t", "z");
    fprintf(fout22, "%s\t\t\t", "z");
	s=LineNumber(fspecies, 10000);
	printf("Species list imported from: %s\n", SPECIES_LIST);
	fclose(fspecies);
	fspecies=fopen(SPECIES_LIST, "r");
	struct Molecule species[s];
	temp=fgets(dataline, 10000, fspecies); /* Read in the header line */
	i=0;
	while (fgets(dataline, 10000, fspecies) != NULL )
	{
		sscanf(dataline, "%s %s %d %d %le %lf %d %lf %le", (species+i)->name, (species+i)->type, &((species+i)->num), &((species+i)->mass), &((species+i)->mix), &((species+i)->upper), &((species+i)->lowertype), &((species+i)->lower), &((species+i)->lower1));
		printf("%8s %2s %4d %4.2d %10.2e %5.2f %3d %5.2f %10.2e\n",(species+i)->name, (species+i)->type, (species+i)->num, (species+i)->mass, (species+i)->mix, (species+i)->upper, (species+i)->lowertype, (species+i)->lower, (species+i)->lower1);
        if (strcmp("X",species[i].type)==0) {numx=numx+1; fprintf(fout21, "%s\t\t\t", (species+i)->name);}
		if (strcmp("F",species[i].type)==0) {numf=numf+1; fprintf(fout22, "%s\t\t\t", (species+i)->name);}
		if (strcmp("C",species[i].type)==0) {numc=numc+1;}
		if (strcmp("A",species[i].type)==0) {numx=numx+1; fprintf(fout21, "%s\t\t\t", (species+i)->name); numa=numa+1;}
		i=i+1;
    }
    fclose(fspecies);
	fprintf(fout21, "%s\n\n", "Air");
	fprintf(fout22, "%s\n\n", "Air");
    fclose(fout21);
    fclose(fout22);
	nums=numx+numf+numc;
	nn=zbin*numx;
	double Con[nn+1], fvec[nn+1], fvec1[nn+1], ConC[zbin*numc+1], Conf[zbin*numf+1]; 
	printf("%s\n", "The species list is imported.");
	printf("%s %d\n", "Number of species in model:", nums);
	printf("%s %d\n", "Number of species to be solved in full:", numx);
	printf("%s %d\n", "In which the number of aerosol species is:", numa);
	printf("%s %d\n", "Number of species to be solved in photochemical equil:", numf);
	printf("%s %d\n", "Number of species assumed to be constant:", numc);
	
	
printf("%s\n\n",fillmi);
	/* Variable initialization */
	int labelx[numx+1], labelc[numc+1], labelf[numf+1], MoleculeM[numx+1], listFix[numx+1], listAER[numa+1], AERCount=1; /* Standard number list of species */
	double Upflux[numx+1], Loflux[numx+1], Depo[numx+1], ConFix[numx+1], mixtemp;
	FILE *fimport;
	FILE *fimportcheck;
	for (i=0; i<s; i++) {
		if (strcmp("X",species[i].type)==0 || strcmp("A",species[i].type)==0) {
			if (IMODE==1) {
				for (j=1; j<=zbin; j++) { Con[(j-1)*numx+numx1]=MM[j]*species[i].mix;} /* Initialized the variables */
			}
			labelx[numx1]=species[i].num;
			if (species[i].num==7) {
				waternum=numx1;
				waterx=1;
			}
			MoleculeM[numx1]=species[i].mass;
			if (species[i].lowertype==1) {
				ConFix[numx1]=species[i].lower1*MM[1];
				Con[numx1]=ConFix[numx1];
				listFix[numx1]=1;
			} else {
				listFix[numx1]=0;
			}
			if (strcmp("A",species[i].type)==0) {
				listAER[AERCount]=numx1;
				AERCount = AERCount+1;
				printf("%s %d\n", "The aerosol species is", numx1);
			}
			numx1=numx1+1;
		}
		if (strcmp("F",species[i].type)==0) {
			labelf[numf1]=species[i].num;
			if (IMODE==1) {
				for (j=1; j<=zbin; j++) { Conf[(j-1)*numf+numf1]=MM[j]*species[i].mix;} /* Initialized the variables */
			}
			numf1=numf1+1;
		}
		if (strcmp("C",species[i].type)==0) {
			labelc[numc1]=species[i].num;
			for (j=1; j<=zbin; j++) {
				ConC[(j-1)*numc+numc1]=MM[j]*species[i].mix; /* Initialize the constants */
			}
			/* import constant mixing ratio list for H2O */
			if (IFIMPORTH2O == 1 && species[i].num == 7) {
				fimport=fopen("Data/ConstantMixing.dat", "r");
				fimportcheck=fopen("Data/ConstantMixingH2O.dat", "w");
				temp=fgets(dataline, 10000, fimport); /* Read in the header line */
				for (j=1; j<=zbin; j++) {
					fscanf(fimport, "%lf\t", &temp1);
					fscanf(fimport, "%le\t", &mixtemp);
					fscanf(fimport, "%le\t", &temp1);
					ConC[(j-1)*numc+numc1]=mixtemp * MM[j];
					fprintf(fimportcheck, "%f\t%e\t%e\n", zl[j], mixtemp, ConC[(j-1)*numc+numc1]);
				}
				fclose(fimport);
				fclose(fimportcheck);
			}
			/* import constant mixing ratio list for H2O */
			if (IFIMPORTCO2 == 1 && species[i].num == 52) {
				fimport=fopen("Data/ConstantMixing.dat", "r");
				fimportcheck=fopen("Data/ConstantMixingCO2.dat", "w");
				temp=fgets(dataline, 10000, fimport); /* Read in the header line */
				for (j=1; j<=zbin; j++) {
					fscanf(fimport, "%lf\t", &temp1);
					fscanf(fimport, "%le\t", &temp1);
					fscanf(fimport, "%le\t", &mixtemp);
					ConC[(j-1)*numc+numc1]=mixtemp * MM[j];
					fprintf(fimportcheck, "%f\t%e\t%e\n", zl[j], mixtemp, ConC[(j-1)*numc+numc1]);
				}
				fclose(fimport);
				fclose(fimportcheck);
			}
			numc1=numc1+1;
		}
    }
	int labels[numx+numf+1];
	
printf("%s\n\n",fillmi);
	
	/* Get Reaction List */
	fzone=fopen(REACTION_LIST, "r");
	s=LineNumber(fzone, 10000);
	fclose(fzone);
	fzone=fopen(REACTION_LIST, "r");
	struct Reaction React[s];
	temp=fgets(dataline, 10000, fzone); /* Read in the header line */
	i=0;
	while (fgets(dataline, 10000, fzone) != NULL )
	{
		sscanf(dataline, "%d %s %d", &((React+i)->dum), (React+i)->type, &((React+i)->num));
		//printf("%d %s %d\n", (React+i)->dum, React[i].type, React[i].num); //ms2022: TMI
		if (strcmp("R",React[i].type)==0) {numr=numr+1;}
		if (strcmp("M",React[i].type)==0) {numm=numm+1;}
		if (strcmp("P",React[i].type)==0) {nump=nump+1;}
		if (strcmp("T",React[i].type)==0) {numt=numt+1;}
		i=i+1;
	}
	fclose(fzone);
	int zone_r[numr+1], zone_m[numm+1], zone_p[nump+1], zone_t[numt+1];
	for (i=0; i<s; i++) {
		if (strcmp("R",React[i].type)==0) {
			zone_r[numr1]=React[i].num;
			numr1=numr1+1;
		}
		;if (strcmp("M",React[i].type)==0) {
			zone_m[numm1]=React[i].num;
			numm1=numm1+1;
		}
		if (strcmp("P",React[i].type)==0) {
			zone_p[nump1]=React[i].num;
			nump1=nump1+1;
		}
		if (strcmp("T",React[i].type)==0) {
			zone_t[numt1]=React[i].num;
			numt1=numt1+1;
		}
	}
	printf("%s\n", "The reaction lists are imported.");
	printf("%s %d\n", "Number of bi-molecular reactions:", numr);
	printf("%s %d\n", "Number of tri-molecular reactions:", numm);
	printf("%s %d\n", "Number of photolysis:", nump);
	printf("%s %d\n", "Number of thermo-dissociations:", numt);
	GetReaction();
	
printf("%s\n\n",fillmi);
	/* get the cross sections and quantum yields of molecules */   
	cross=dmatrix(1,nump,0,NLAMBDA-1);
	crosst=dmatrix(1,nump,0,NLAMBDA-1);
	qy=dmatrix(1,nump,0,NLAMBDA-1);
	qyt=dmatrix(1,nump,0,NLAMBDA-1);
    char dirroutep[1024]="Library/PhotochemOpa/";
    char photochemfile[1024];
//ms23	int stdcross[nump+1];
        stdcross=ivector(0,nump); //try this instead
//ms23	double qysum[nump+1];
        qysum=dvector(0,nump); //try this instead
	fcheck=fopen("AuxillaryOut/CrossSectionCheck.dat","w");
//        printf("%s\n","ReactionP: "); //ms2021: check
	for (i=1; i<=nump; i++) {
//            printf("%.2f ",ReactionP[zone_p[i]][8]); //ms2021: check
		stdcross[i]=ReactionP[zone_p[i]][1];
		qytype=ReactionP[zone_p[i]][8];
		qysum[i]=ReactionP[zone_p[i]][7];
		j=0;
		while (species[j].num != stdcross[i]) {j=j+1;}
		/* printf("%s\n",species[j].name); */
        strcpy(photochemfile,dirroutep);
        strcat(photochemfile,species[j].name);
		fp=fopen(photochemfile, "r");
		fp1=fopen(photochemfile, "r");
		s=LineNumber(fp, 1000);
		/* printf("%d\n",s); */
		wavep=dvector(0,s-1);
		crossp=dvector(0,s-1);
		qyp=dvector(0,s-1);
		qyp1=dvector(0,s-1);
		qyp2=dvector(0,s-1);
		qyp3=dvector(0,s-1);
		qyp4=dvector(0,s-1);
		qyp5=dvector(0,s-1);
		qyp6=dvector(0,s-1);
		qyp7=dvector(0,s-1);
		crosspt=dvector(0,s-1);
		qypt=dvector(0,s-1);
		qyp1t=dvector(0,s-1);
		qyp2t=dvector(0,s-1);
		qyp3t=dvector(0,s-1);
		qyp4t=dvector(0,s-1);
		qyp5t=dvector(0,s-1);
		qyp6t=dvector(0,s-1);
		qyp7t=dvector(0,s-1);
        k=0;
        if (qytype==1) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf", wavep+k, crossp+k, crosspt+k, qyp+k, qypt+k);
				k=k+1; }
        }
        if (qytype==2) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp+k, qypt+k);
				k=k+1; }
        }
        if (qytype==3) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp+k, qypt+k);
				k=k+1; }
        }
	if (qytype==4) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp+k, qypt+k);
				k=k+1; }
        }
	if (qytype==5) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp+k, qypt+k);
				k=k+1; }
        }
	if (qytype==6) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp+k, qypt+k);
				k=k+1; }
        }
	if (qytype==7) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp6+k, qyp6t+k, qyp+k, qypt+k);
				k=k+1; }
        }
	if (qytype==8) {
			while (fgets(dataline, 1000, fp1) != NULL ) {
				sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp6+k, qyp6t+k, qyp7+k, qyp7t+k, qyp+k, qypt+k);
				k=k+1; }
        }
		fclose(fp);
		fclose(fp1);
		Interpolation(wavelength, NLAMBDA, *(cross+i), wavep, crossp, s, 0);
		Interpolation(wavelength, NLAMBDA, *(qy+i), wavep, qyp, s, 0);
		Interpolation(wavelength, NLAMBDA, *(crosst+i), wavep, crosspt, s, 0);
		Interpolation(wavelength, NLAMBDA, *(qyt+i), wavep, qypt, s, 0);
		free_dvector(wavep,0,s-1);
		free_dvector(crossp,0,s-1);
		free_dvector(qyp,0,s-1);
		free_dvector(qyp1,0,s-1);
                free_dvector(qyp2,0,s-1);
		free_dvector(qyp3,0,s-1);
		free_dvector(qyp4,0,s-1);
		free_dvector(qyp5,0,s-1);
		free_dvector(qyp6,0,s-1);
		free_dvector(qyp7,0,s-1);
		free_dvector(crosspt,0,s-1);
		free_dvector(qypt,0,s-1);
		free_dvector(qyp1t,0,s-1);
                free_dvector(qyp2t,0,s-1);
		free_dvector(qyp3t,0,s-1);
		free_dvector(qyp4t,0,s-1);
		free_dvector(qyp5t,0,s-1);
		free_dvector(qyp6t,0,s-1);
		free_dvector(qyp7t,0,s-1);
		printf("%s %s %s\n", "The", species[j].name, "Cross section and quantum yield data are imported.");
		fprintf(fcheck, "%s %s %s\n", "The", species[j].name, "Cross section and quantum yield data are imported.");
		for (j=0;j<NLAMBDA;j++) {fprintf(fcheck, "%lf %le %le %lf %lf\n", wavelength[j], cross[i][j], crosst[i][j], qy[i][j], qyt[i][j]);} 
	}
	
	/* cross section of aerosols */
	double *crossp1, *crossp2, *crossp3;
	double crossw1[NLAMBDA], crossw2[NLAMBDA], crossw3[NLAMBDA];
	fp=fopen(AERRADFILE1,"r");
	fp1=fopen(AERRADFILE1,"r");
	s=LineNumber(fp, 1000);
	wavep=dvector(0,s-1);
	crossp1=dvector(0,s-1);
	crossp2=dvector(0,s-1);
	crossp3=dvector(0,s-1);
	k=0;
	while (fgets(dataline, 1000, fp1) != NULL ) {
		sscanf(dataline, "%lf %lf %lf %lf", wavep+k, crossp1+k, crossp2+k, crossp3+k);
		k=k+1; 
	}
	fclose(fp);
	fclose(fp1);
	Interpolation(wavelength, NLAMBDA, crossw1, wavep, crossp1, s, 0);
	Interpolation(wavelength, NLAMBDA, crossw2, wavep, crossp2, s, 0);
	Interpolation(wavelength, NLAMBDA, crossw3, wavep, crossp3, s, 0);
	free_dvector(wavep,0,s-1);
	free_dvector(crossp1,0,s-1);
	free_dvector(crossp2,0,s-1);
	free_dvector(crossp3,0,s-1);
	for (i=0; i<NLAMBDA; i++) {
		crossa[1][i] = crossw1[i];
		sinab[1][i]  = crossw2[i]/(crossw1[i]+1.0e-24);
		asym[1][i]   = crossw3[i];
	}
	fp=fopen(AERRADFILE2,"r");
	fp1=fopen(AERRADFILE2,"r");
	s=LineNumber(fp, 1000);
	wavep=dvector(0,s-1);
	crossp1=dvector(0,s-1);
	crossp2=dvector(0,s-1);
	crossp3=dvector(0,s-1);
	k=0;
	while (fgets(dataline, 1000, fp1) != NULL ) {
		sscanf(dataline, "%lf %lf %lf %lf", wavep+k, crossp1+k, crossp2+k, crossp3+k);
		k=k+1; 
	}
	fclose(fp);
	fclose(fp1);
	Interpolation(wavelength, NLAMBDA, crossw1, wavep, crossp1, s, 0);
	Interpolation(wavelength, NLAMBDA, crossw2, wavep, crossp2, s, 0);
	Interpolation(wavelength, NLAMBDA, crossw3, wavep, crossp3, s, 0);
	free_dvector(wavep,0,s-1);
	free_dvector(crossp1,0,s-1);
	free_dvector(crossp2,0,s-1);
	free_dvector(crossp3,0,s-1);
	for (i=0; i<NLAMBDA; i++) {
		crossa[2][i] = crossw1[i];
		sinab[2][i]  = crossw2[i]/(crossw1[i]+1.0e-24);
		asym[2][i]   = crossw3[i];
	}
	printf("%s\n", "Cross sections of the aerosol are imported.");
	fprintf(fcheck, "%s\n", "Cross sections of the aerosol are imported.");
	for (j=0;j<NLAMBDA;j++) {fprintf(fcheck, "%lf %e %e %f %f %f %f\n", wavelength[j], crossa[1][j], crossa[2][j], sinab[1][j], sinab[2][j], asym[1][j], asym[2][j]);}
	fclose(fcheck);
	
    double importpres, importtemp, importmm; 
    double importpres_save[zbin], importnn, importnn_save[zbin+1][NSP+1], importnn_save1[zbin], importnn_save2[zbin], importpl[zbin]; //rh2024, change to interpolation when importing
    
    char dirroute[1024];
    char atomfile[1024];
    char outstdt[1024];
    char outnewtemp[1024];
    char crossfile[1024];
    
    char outchemstatus[1024];
    char outoldtemp[1024];
    char outstd0[1024];
    char outstd[1024];
    char outphotorate[1024];
    char outfile1[1024];
    char outfile2[1024];
    char outchemicalrate[1024];
    char outcolumn[1024];
    char outradiation[1024];
    char outconvergence[1024];
    char outtimescale[1024];
    char outbalance[1024];
    char outhistory[1024];
    char outoldstd[1024];

    char outrtdiag[1024]; //for RadTrans diagnostics
    char outrcdiag[1024]; //for Rad-Conv diagnostics
    char outcondiag[1024]; //for Condensibles diagnostics
    
    strcpy(dirroute,OUT_DIR);
    strcpy(atomfile,dirroute);
    strcat(atomfile,"atom.dat");
    printf("%s\t%s\n","Prepare to get atom abundance from", atomfile);
    
    double totalmass;
    double totalnumber;
    double totalmix;
    double heliumnumber;
    double TVARTOTAL;
    
        //for (j=1; j<=zbin; j++)  for (i=1; i<=NSP; i++) clouds[j][i]=0.0; //initializing for Climate

        /* IMODE = 4: from existing files */
    if (IMODE == 4) {
        strcpy(outstd0,dirroute);
        strcat(outstd0,"/ConcentrationSTD_C.dat");
        printf("%s\t%s\n","Prepare to get initial molecular concentration from", outstd0);
        fimport=fopen(outstd0, "r");
        fimportcheck=fopen("AuxillaryOut/Fimportcheck.dat","w");
        temp=fgets(dataline, 10000, fimport); // Read in the header line //
        temp=fgets(dataline, 10000, fimport); // Read in the header line //
        for (j=1; j<=zbin; j++) {
            fscanf(fimport, "%lf\t", &temp1);
            fprintf(fimportcheck, "%lf\t", temp1);
            fscanf(fimport, "%lf\t", &temp1);
            fscanf(fimport, "%lf\t", &temp1);
            fscanf(fimport, "%lf\t", &importtemp);
            fscanf(fimport, "%le\t", &importpres);
            // printf("%f %e\n",importtemp,importpres);
            importpres_save[zbin-j]=log(importpres);
            importmm = importpres/KBOLTZMANN/importtemp*1.0E-6;
            //rh2024 change to interpolate with respect to pres
            for (i=1; i<=NSP; i++) {
                // fscanf(fimport, "%le\t", &xx[j][i]);
                // xx[j][i]  = xx[j][i]/importmm*MM[j];
                fscanf(fimport, "%le\t", &importnn);
                importnn_save[j][i]=log(importnn);
                fprintf(fimportcheck, "%e\t", exp(importnn_save[j][i]));
            }
            fprintf(fimportcheck,"\n");
        }
        fclose(fimport);
        fclose(fimportcheck);
        //rh2024 change to interpolate with respect to pres
        for (i=1; i<=NSP; i++) {
            for (j=1; j<=zbin; j++) {
                importnn_save1[zbin-j] = importnn_save[j][i];
                importpl[zbin-j] = log(pl[j]);
            }
            Interpolation(importpl,zbin,importnn_save2,importpres_save,importnn_save1,zbin,0);
            if (importnn_save2[0]==0) { importnn_save2[0]=importnn_save1[0]-importpres_save[0]+importpl[0];
            }
            if (importnn_save2[zbin-1]==0) { importnn_save2[zbin-1]=importnn_save1[zbin-1]-importpres_save[zbin-1]+importpl[zbin-1];
            }
            if (i==43) {
                for (j=1; j<=zbin; j++) {
                    printf("%f %f %f %f\n",importpl[j-1],importnn_save2[j-1],importpres_save[j-1],importnn_save1[j-1]);
                }
            }
            for (j=1; j<=zbin; j++) {
                xx[j][i]=exp(importnn_save2[zbin-j]);
                if isnan(xx[j][i]) xx[j][i]=0.0;
            }
        }
        Convert2(Con, ConC, Conf, labelx, labelc, labelf); //getting Con from XX
    
    //strcpy(outstd,dirroute);
    //strcat(outstd,"/ConcentrationSTD_T.dat");
    //printf("%s\t%s\n","Prepare to write molecular concentration to", outstd);
    printf("%s\n", "Variable initialization successful");
    //printout_std(z,outstd);
    
    strcpy(outoldstd,dirroute);
    strcat(outoldstd,"/ConcentrationSTD_Old.dat");
    printout_std(z,outoldstd);
    }

    if (IMODE < 4) {

    /* compute the initial molecular abundances */
    for (i=1; i<=numx; i++) {labels[i]=labelx[i];}
    for (i=1; i<=numf; i++) {labels[numx+i]=labelf[i];}
    mixequil=dmatrix(1,zbin,1,numx+numf);
    for (j=1; j<=zbin; j++) {
        for (i=1; i<=numx+numf; i++) {
            mixequil[j][i]=0.0;
        }
    }
    chemquil(pl, tl, zbin+1, labels, numx+numf, mixequil,atomfile);
    checkmixequil(numx+numf, mixequil);
    for (j=1; j<=zbin; j++) {
        for (i=1; i<=numx; i++) {
            Con[(j-1)*numx+i]=MM[j]*mixequil[j][i];
            //printf("X %d %d %e\n", j, labelx[i], mixequil[j][i]); //ms2022: TMI
        }
        for (i=1; i<=numf; i++) {
            Conf[(j-1)*numf+i]=MM[j]*mixequil[j][i+numx];
            //printf("F %d %d %e\n", j, labelf[i], mixequil[j][i+numx]); //ms2022: TMI
        }
    }
    free_dmatrix(mixequil,1,zbin,1,numx+numf);
    
    /* Generate General Variables */
    Convert1(Con, ConC, Conf, labelx, labelc, labelf); //getting XX from Con
    printf("%s\n", "Variable initialization successful");
    } //IMODE<4

    /* Generate Mean Molecular Mass */
    for (j=1; j<=zbin; j++) {
        totalnumber=0.0;
        totalmass=0.0;
        for (i=1; i<=numx; i++) {
            totalnumber += Con[(j-1)*numx+i];
            totalmass   += Con[(j-1)*numx+i]*MoleculeM[i];
        }
        /* count for Helium */
        heliumnumber = fmax(MM[j] - totalnumber,0.0);
        totalnumber += heliumnumber;
        totalmass += heliumnumber*4.0;
        meanmolecular[j] = totalmass/totalnumber;
    }
    printf("%s %d %s %2.2f\n","Mean Molecular Mass at layer ", j-1, "is", meanmolecular[j-1]);
    printf("%s %d %s %2.2e\n","Helium mixing ratio at layer ", j-1, "is", heliumnumber/totalnumber);
    
    /* Obtain the opacity */
    opacCO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2O = dmatrix(1,zbin,0,NLAMBDA-1);
    opacOH = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2CO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2O2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacO3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCH4 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H4 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H6 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCH2O2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHCN = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNH3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacN2O = dmatrix(1,zbin,0,NLAMBDA-1);
    opacN2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHNO3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2S = dmatrix(1,zbin,0,NLAMBDA-1);
    opacSO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacOCS = dmatrix(1,zbin,0,NLAMBDA-1);
    
    readcia();
    planckmeanCIA();
    printf("CIA mean opacity in the infrared calculated!\n");
    
//clock_gettime(CLOCK_REALTIME, &tstart);
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacC2H2.dat");
    readcross(crossfile, opacC2H2);
    planckmean(MeanC2H2, SMeanC2H2, opacC2H2);
    
//clock_gettime(CLOCK_REALTIME, &tend);
//t_passed = ((double)tend.tv_sec*1e9 + tend.tv_nsec) - ((double)tstart.tv_sec*1e9 + tstart.tv_nsec);
//printf("\n%s %f\n","reading ONE opacity file took [ms]:", t_passed/1e6); 

    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacC2H4.dat");
    readcross(crossfile, opacC2H4);
    planckmean(MeanC2H4, SMeanC2H4, opacC2H4);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacC2H6.dat");
    readcross(crossfile, opacC2H6);
    planckmean(MeanC2H6, SMeanC2H6, opacC2H6);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacCH2O2.dat");
    readcross(crossfile, opacCH2O2);
    planckmean(MeanCH2O2, SMeanCH2O2, opacCH2O2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacCH4.dat");
    readcross(crossfile, opacCH4);
    planckmean(MeanCH4, SMeanCH4, opacCH4);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacCO2.dat");
    readcross(crossfile, opacCO2);
    planckmean(MeanCO2, SMeanCO2, opacCO2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacCO.dat");
    readcross(crossfile, opacCO);
    planckmean(MeanCO, SMeanCO, opacCO);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacH2CO.dat");
    readcross(crossfile, opacH2CO);
    planckmean(MeanH2CO, SMeanH2CO, opacH2CO);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacH2O2.dat");
    readcross(crossfile, opacH2O2);
    planckmean(MeanH2O2, SMeanH2O2, opacH2O2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacH2O.dat");
    readcross(crossfile, opacH2O);
    planckmean(MeanH2O, SMeanH2O, opacH2O);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacH2S.dat");
    readcross(crossfile, opacH2S);
    planckmean(MeanH2S, SMeanH2S, opacH2S);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacHCN.dat");
    readcross(crossfile, opacHCN);
    planckmean(MeanHCN, SMeanHCN, opacHCN);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacHNO3.dat");
    readcross(crossfile, opacHNO3);
    planckmean(MeanHNO3, SMeanHNO3, opacHNO3);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacHO2.dat");
    readcross(crossfile, opacHO2);
    planckmean(MeanHO2, SMeanHO2, opacHO2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacN2O.dat");
    readcross(crossfile, opacNO2);
    planckmean(MeanNO2, SMeanNO2, opacNO2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacN2.dat");
    readcross(crossfile, opacN2);
    planckmean(MeanN2, SMeanN2, opacN2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacNH3.dat");
    readcross(crossfile, opacNH3);
    planckmean(MeanNH3, SMeanNH3, opacNH3);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacNO2.dat");
    readcross(crossfile, opacNO2);
    planckmean(MeanNO2, SMeanNO2, opacNO2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacNO.dat");
    readcross(crossfile, opacNO);
    planckmean(MeanNO, SMeanNO, opacNO);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacO2.dat");
    readcross(crossfile, opacO2);
    planckmean(MeanO2, SMeanO2, opacO2);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacO3.dat");
    readcross(crossfile, opacO3);
    planckmean(MeanO3, SMeanO3, opacO3);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacOCS.dat");
    readcross(crossfile, opacOCS);
    planckmean(MeanOCS, SMeanOCS, opacOCS);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacOH.dat");
    readcross(crossfile, opacOH);
    planckmean(MeanOH, SMeanOH, opacOH);
    
    strcpy(crossfile,CROSSHEADING);
    strcat(crossfile,"opacSO2.dat");
    readcross(crossfile, opacSO2);
    planckmean(MeanSO2, SMeanSO2, opacSO2);
    
//clock_gettime(CLOCK_REALTIME, &tend);
//t2_passed = ((double)tend.tv_sec*1e9 + tend.tv_nsec) - ((double)tstart.tv_sec*1e9 + tstart.tv_nsec);
//printf("\n%s %f\n","reading ALL opacity files took [ms]:", t2_passed/1e6); 

    /* Print Out Initialization */
    strcpy(outstdt,dirroute);
    strcat(outstdt,"/ConcentrationSTD_T.dat");
    printf("%s\t%s\n","Prepare to print out thermal STD file to", outstdt);
    printout_std_t(z,outstdt);
    
    /* New Temperature - Initialization */
    strcpy(outnewtemp,dirroute);
    strcat(outnewtemp,"/NewTemperature.dat"); //ms2021: "/" removed
    printf("%s\t%s\n","Prepare to print out new temperature to", outnewtemp);
    
    /* New Temperature - Initialization */
    strcpy(outrtdiag,dirroute);
    strcat(outrtdiag,"/Diagnostic_RT.dat"); //ms2021: "/" removed
    printf("%s\t%s\n","Prepare to print out RT diagnostics to", outrtdiag);
    
    /* New Temperature - Initialization */
    strcpy(outrcdiag,dirroute);
    strcat(outrcdiag,"/Diagnostic_RC.dat"); //ms2021: "/" removed
    printf("%s\t%s\n","Prepare to print out RC diagnostics to", outrcdiag);
    
    /* New Temperature - Initialization */
    strcpy(outcondiag,dirroute);
    strcat(outcondiag,"/Diagnostic_condens.dat"); //ms2021: "/" removed
    printf("%s\t%s\n","Prepare to print out condensibles diagnostics to", outcondiag);
    
//exit(0);//ms: for testing purposes

    if(RadConv_Solver == 0) GreyTemp(P,outnewtemp,TINTSET); 
    
    if(RadConv_Solver == 1) {
        ms_Climate(tempeq,P,T,TINTSET,outnewtemp,outrtdiag,outrcdiag,outcondiag); //ms22: 2stream switch added
        for (j=0; j<=zbin; j++) Tnew[j]=tempeq[j];
    }
//**************************************************************
    /* Iteration */
//**************************************************************
    for (i=0; i<NMAX; i++) {
        
        /* Determine if converged */
        TVARTOTAL = 0.0;
        for (j=0; j<=zbin; j++) {
            TVARTOTAL += fabs(Tnew[j] - T[j]);
        }
        TVARTOTAL /= zbin;
        printf("%s %f %s\n", "Temperature variation is ", TVARTOTAL, "K");
        
        /* Check convergence */
        if (TVARTOTAL<1.0) {
            printf("%s\n", "EPACRIS converged!");
            /* fprintf(fstat,"%s\n", "converged!");*/
            printout_std_t(z,outstdt);
            break;
        }
        
        /* If not converged, update */
        /* update */
        for (j=0; j<=zbin; j++) {
            T[j] = Tnew[j];
            Tdoub[2*j] = T[j]; //ms2023: double grid
        }
        for (j=1; j<=zbin; j++) {
            tl[j] = (T[j]+T[j-1])/2.0;
            Tdoub[2*j-1] = tl[j]; //ms2023: double grid
        }
        z[0] = 0.0;
        for (j=1; j<=zbin; j++) {
            scaleheight = KBOLTZMANN * tl[j] / meanmolecular[j] / AMU / GA /1000.0 ; /* km */
            z[j] = z[j-1] - scaleheight*log(P[j]/P[j-1]);
            zdoub[2*j] = z[j]; //ms2023: double grid
        }
        for (j=1; j<=zbin; j++) {
            zl[j] = (z[j]+z[j-1])/2.0;
            zdoub[2*j-1] = zl[j]; //ms2023: double grid
        }
        for (j=0; j<=zbin; j++) {
            MMZ[j] = P[j]/KBOLTZMANN/T[j]*1.0E-6;
            MMdoub[2*j] = MMZ[j]; //ms2023: double grid
        }
        for (j=1; j<=zbin; j++) {
            MM[j]=pl[j]/KBOLTZMANN/tl[j]*1.0E-6; /*unit: Molecule cm-3*/
            MMdoub[2*j-1] = MM[j]; //ms2023: double grid
            printf("%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
        }
        /* Print out */
        printout_std_t(z,outstdt);

        printf("%s\n", "The Z-T-P profile is re-calculated.");
        
        /* Re-calculate the initial mixing ratio from chemical equilibrium  */
        //if (IMODE==0 || IMODE==4) {
        if (IMODE==0) {
            mixequil=dmatrix(1,zbin,1,numx+numf);
            for (j=1; j<=zbin; j++) {
                for (ii=1; ii<=numx+numf; ii++) {
                    mixequil[j][ii]=0.0;
                }
            }
            chemquil(pl, tl, zbin+1, labels, numx+numf, mixequil, atomfile);
            checkmixequil(numx+numf, mixequil);
            for (j=1; j<=zbin; j++) {
                totalmix = 0.0;
                for (ii=1; ii<=numx+numf; ii++) {
                    xx[j][labels[ii]]=MM[j]*mixequil[j][ii];
                    totalmix += mixequil[j][ii];
                }
                printf("%s %d %s %2.2f\n","Total mixing ratio at layer", j, "is", totalmix);
            }
            free_dmatrix(mixequil,1,zbin,1,numx+numf);
        Convert2(Con, ConC, Conf, labelx, labelc, labelf);
        printf("%s\n", "The thermochem equilibrium composition is re-calculated.");
        } // Imode==0
        
        if (IMODE==4) {
        Convert2(Con, ConC, Conf, labelx, labelc, labelf);
        printf("%s\n", "Con[i] re-calculated to incorporate any potential XX[i] changes from previous climate calculations");
        }
        
        /* Update Mean Molecular Mass */
        for (j=1; j<=zbin; j++) {
            totalnumber=0.0;
            totalmass=0.0;
            for (ii=1; ii<=numx; ii++) {
                totalnumber += Con[(j-1)*numx+ii];
                totalmass   += Con[(j-1)*numx+ii]*MoleculeM[ii];
            }
            /* count for Helium */
            heliumnumber = fmax(MM[j] - totalnumber,0.0);
            totalnumber += heliumnumber;
            totalmass += heliumnumber*4.0;
            meanmolecular[j] = totalmass/totalnumber;
            printf("%s %d %s %2.2f\n","Mean Molecular Mass at layer", j, "is", meanmolecular[j]);
            printf("%s %d %s %2.2e\n","Helium mixing ratio at layer", j, "is", heliumnumber/totalnumber);
        }
        
        /* Update Opacities */
        readcia();
        planckmeanCIA();
        printf("CIA mean opacity in the infrared calculated!\n");
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacCO2.dat");
        readcross(crossfile, opacCO2);
        planckmean(MeanCO2, SMeanCO2, opacCO2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacO2.dat");
        readcross(crossfile, opacO2);
        planckmean(MeanO2, SMeanO2, opacO2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacSO2.dat");
        readcross(crossfile, opacSO2);
        planckmean(MeanSO2, SMeanSO2, opacSO2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacH2O.dat");
        readcross(crossfile, opacH2O);
        planckmean(MeanH2O, SMeanH2O, opacH2O);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacOH.dat");
        readcross(crossfile, opacOH);
        planckmean(MeanOH, SMeanOH, opacOH);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacH2CO.dat");
        readcross(crossfile, opacH2CO);
        planckmean(MeanH2CO, SMeanH2CO, opacH2CO);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacH2O2.dat");
        readcross(crossfile, opacH2O2);
        planckmean(MeanH2O2, SMeanH2O2, opacH2O2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacHO2.dat");
        readcross(crossfile, opacHO2);
        planckmean(MeanHO2, SMeanHO2, opacHO2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacH2S.dat");
        readcross(crossfile, opacH2S);
        planckmean(MeanH2S, SMeanH2S, opacH2S);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacCO.dat");
        readcross(crossfile, opacCO);
        planckmean(MeanCO, SMeanCO, opacCO);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacO3.dat");
        readcross(crossfile, opacO3);
        planckmean(MeanO3, SMeanO3, opacO3);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacCH4.dat");
        readcross(crossfile, opacCH4);
        planckmean(MeanCH4, SMeanCH4, opacCH4);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacNH3.dat");
        readcross(crossfile, opacNH3);
        planckmean(MeanNH3, SMeanNH3, opacNH3);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacC2H2.dat");
        readcross(crossfile, opacC2H2);
        planckmean(MeanC2H2, SMeanC2H2, opacC2H2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacC2H4.dat");
        readcross(crossfile, opacC2H4);
        planckmean(MeanC2H4, SMeanC2H4, opacC2H4);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacC2H6.dat");
        readcross(crossfile, opacC2H6);
        planckmean(MeanC2H6, SMeanC2H6, opacC2H6);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacCH2O2.dat");
        readcross(crossfile, opacCH2O2);
        planckmean(MeanCH2O2, SMeanCH2O2, opacCH2O2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacHCN.dat");
        readcross(crossfile, opacHCN);
        planckmean(MeanHCN, SMeanHCN, opacHCN);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacNO.dat");
        readcross(crossfile, opacNO);
        planckmean(MeanNO, SMeanNO, opacNO);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacNO2.dat");
        readcross(crossfile, opacNO2);
        planckmean(MeanNO2, SMeanNO2, opacNO2);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacHNO3.dat");
        readcross(crossfile, opacHNO3);
        planckmean(MeanHNO3, SMeanHNO3, opacHNO3);
        
        strcpy(crossfile,CROSSHEADING);
        strcat(crossfile,"opacOCS.dat");
        readcross(crossfile, opacOCS);
        planckmean(MeanOCS, SMeanOCS, opacOCS);
        
        /* New Temperature */
        if(RadConv_Solver == 0) GreyTemp(P,outnewtemp,TINTSET);
        
        if(RadConv_Solver == 1) {
            ms_Climate(tempeq,P,T,TINTSET,outnewtemp,outrtdiag,outrcdiag,outcondiag);
            for (j=0; j<=zbin; j++) Tnew[j]=tempeq[j];
        }
        
    }
    
//**************************************************************
    /* fprintf(fstat,"%s %f\n", "Temperature variation is", TVARTOTAL);
    fclose(fstat); */
    
    /* Clean up */
    free_dmatrix(opacCO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2O,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacOH,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2CO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2O2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacO3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCH4,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H4,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H6,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCH2O2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHCN,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNH3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacN2O,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHNO3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2S,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacSO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacOCS,1,zbin,0,NLAMBDA-1);

    /* Clean up */
    free_dmatrix(cross,1,nump,0,NLAMBDA-1);
    free_dmatrix(qy,1,nump,0,NLAMBDA-1);
    free_dmatrix(crosst,1,nump,0,NLAMBDA-1);
    free_dmatrix(qyt,1,nump,0,NLAMBDA-1);

	
//=== END: MAIN PROGRAM ==================================
//TRANSFER gattaca output into folder:
//printf("%s\n","IS THIS EXECUTED?");
if (argc > 1)
{
    //printf("%s\n","AND THIS?");
    char cmd[500];
    char * jobid = strtok(argv[1],".");
    sprintf(cmd,"%s %s%s%s %s%s","mv",argv[2],".o",jobid,OUT_DIR,"/" );
    system(cmd);
    //now archive input file
    sprintf(cmd,"%s%s %s%s%s%s%s","cp Input/",IN_FILE_NAME,OUT_DIR,"/",jobid,"_",IN_FILE_NAME );
    system(cmd);

    //printf("%s\n","OR THAT?");
}
//========================================================
}
//========================================================
//========================================================

