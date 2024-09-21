#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "planet_chem.h"

#include "constant.h"
#include "routine.h"
#include "global_chem.h"
#include "GetData.c"
#include "chemequil.c"
#include "Interpolation.c"
#include "nrutil.h"
#include "nrutil.c"
#include "RadTransfer.c"
#include "Photorate.c"
#include "ChemReaction.c"
#include "lubksb.c"
#include "ludcmp_mod.c"
#include "Convert.c"
#include "BTridiagonal_mod.c"
#include "TPPara.c"
#include "TPScale.c"
#include "printout_exp.c"
#include "printout_std_t_exp.c"
#include "printoutrate_exp.c"
#include "RefIdx.c"
#include "fallvelocity.c"
#include "saturationpressure.c"
#include "GlobalBalance.c"
#include "GetGibbsForm.c"
#include "ReverseReactionRate.c"
#include "readcia.c"
#include "readcross.c"
#include "planckmean.c"
#include "planckmeanCIA.c"

#include "ReactionRate_Moses.c"
#include "ReactionRateM_Moses.c"
#include "ReactionRateT_Moses.c"
#include "GreyTemp.c"


/* external (global) variables */

double meanmolecular[zbin+1];
double thickl;
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
double **opacC2H2, **opacC2H4, **opacC2H6, **opacCH2O2;
double **opacHCN, **opacN2O, **opacNO, **opacNO2, **opacHNO3, **opacOCS;

double MeanCO2[zbin+1], MeanO2[zbin+1], MeanSO2[zbin+1], MeanH2O[zbin+1], MeanOH[zbin+1], MeanH2CO[zbin+1];
double MeanH2O2[zbin+1], MeanHO2[zbin+1], MeanH2S[zbin+1], MeanCO[zbin+1], MeanO3[zbin+1], MeanCH4[zbin+1];
double MeanNH3[zbin+1];	
double MeanC2H2[zbin+1], MeanC2H4[zbin+1], MeanC2H6[zbin+1], MeanCH2O2[zbin+1];
double MeanHCN[zbin+1], MeanN2O[zbin+1], MeanNO[zbin+1], MeanNO2[zbin+1], MeanOCS[zbin+1], MeanHNO3[zbin+1];

double SMeanCO2[zbin+1], SMeanO2[zbin+1], SMeanSO2[zbin+1], SMeanH2O[zbin+1], SMeanOH[zbin+1], SMeanH2CO[zbin+1];
double SMeanH2O2[zbin+1], SMeanHO2[zbin+1], SMeanH2S[zbin+1], SMeanCO[zbin+1], SMeanO3[zbin+1], SMeanCH4[zbin+1];
double SMeanNH3[zbin+1];
double SMeanC2H2[zbin+1], SMeanC2H4[zbin+1], SMeanC2H6[zbin+1], SMeanCH2O2[zbin+1];
double SMeanHCN[zbin+1], SMeanN2O[zbin+1], SMeanNO[zbin+1], SMeanNO2[zbin+1], SMeanOCS[zbin+1], SMeanHNO3[zbin+1];

double H2H2CIA[zbin+1][NLAMBDA], H2HeCIA[zbin+1][NLAMBDA], H2HCIA[zbin+1][NLAMBDA], N2H2CIA[zbin+1][NLAMBDA], N2N2CIA[zbin+1][NLAMBDA], CO2CO2CIA[zbin+1][NLAMBDA];
double MeanH2H2CIA[zbin+1], MeanH2HeCIA[zbin+1], MeanH2HCIA[zbin+1], MeanN2H2CIA[zbin+1], MeanN2N2CIA[zbin+1], MeanCO2CO2CIA[zbin+1];
double SMeanH2H2CIA[zbin+1], SMeanH2HeCIA[zbin+1], SMeanH2HCIA[zbin+1], SMeanN2H2CIA[zbin+1], SMeanN2N2CIA[zbin+1], SMeanCO2CO2CIA[zbin+1];

double rainoutrate[zbin+1][NSP+1];
double Vesc[NSP+1], VFall[zbin+1];
double nsH2O[zbin+1], nsH2SO4[zbin+1], nsS8[zbin+1], tcondfH2O[zbin+1], tcondfH2SO4[zbin+1], tcondfS8[zbin+1];

double kk[zbin+1][NKin+1], kkM[zbin+1][NKinM+1], kkT[zbin+1][NKinT+1];
double Rkk[zbin+1][NKin+1], RkkM[zbin+1][NKinM+1], RkkT[zbin+1][NKinT+1];
int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
double **DM, **dl, KE[zbin+1];

int    numr=0, numm=0, numt=0, nump=0, numx=0, numc=0, numf=0, numa=0, waternum=0, waterx=0;

double xx[zbin+1][NSP+1];

double mkv[zbin+1], Tnew[zbin+1], Pnew[zbin+1];

double GibbsForm[NSP+1][zbin+1];

main()
{
	  int s,i,ii,j,jj,jjj,k,nn,qytype,stdnum,iradmax,singularflag;
	  int numr1=1,numm1=1,nump1=1,numt1=1;
	  int nums;
	  int numx1=1, numf1=1, numc1=1;
	  char *temp;
	  char dataline[10000];
	  double temp1, wavetemp, crosstemp, DD, GA, DenZ;
	  double z[zbin+1], T[zbin+1], PP[zbin+1], P[zbin+1], zmin, zmax;
      double **JJ, **cross, **qy, *wavep, *crossp, *qyp, *qyp1, *qyp2, *qyp3, *qyp4, *qyp5, *qyp6, *qyp7;
	  double **crosst, **qyt, *crosspt, *qypt, *qyp1t, *qyp2t, *qyp3t, *qyp4t, *qyp5t, *qyp6t, *qyp7t;
	  FILE *fspecies, *fzone, *fhenry, *fp, *fp1, *fp2, *fp3;
      FILE *fout, *fout1, *fout21, *fout22, *fout3, *fout4, *fcheck, *ftemp, *fout5, *foutp, *fcheckgibbs;

      GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET; /* Planet Surface Gravity Acceleration, in SI */	
	
    
    double dlambda, start, interval;
	start = log10(LAMBDALOW);
	interval = log10(LAMBDAHIGH) - log10(LAMBDALOW);
	dlambda = interval / (NLAMBDA-1.0);
	for (i=0; i<NLAMBDA; i++){
		wavelength[i] = pow(10.0, start+i*dlambda)*1.0E+9; /* in nm */
	}
    
      /* Set the wavelength for computation of photolysis rate */
      double refidx0, DenS;
	  iradmax = 0;
	  DenS=101325.0/KBOLTZMANN/273.0*1.0E-6; 
      for (i=0; i<NLAMBDA; i++) {
		  if (wavelength[i]<=WaveMax) { iradmax = i; }
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
	printf("%s %d\n", "The maximum index for the UV-Visible Rad Transfer and Photolysis is", iradmax);
      
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
		if (IFUVMULT == 1) {
			if (wavelength[i] < 200.0) {
				solar[i] = solar[i] * FUVMULT;
			} else if (wavelength[i] < 300.0) {
				solar[i] = solar[i] * MUVMULT;
			} else if (wavelength[i] < 400.0) {
				solar[i] = solar[i] * NUVMULT;
			}
		}
    }
    i=0;
    while (solar[i]>0 || wavelength[i]<9990 ) { i++;}
    for (j=i; j<NLAMBDA; j++) {
        solar[j] = solar[i-1]*pow(wavelength[i-1],4)/pow(wavelength[j],4);
    }
	printf("%s\n", "The stellar radiation data are imported.");
	
	/* Get the species list */
	fspecies=fopen(SPECIES_LIST, "r");
	s=LineNumber(fspecies, 10000);
	printf("Species list: \n");
	fclose(fspecies);
	fspecies=fopen(SPECIES_LIST, "r");
	struct Molecule species[s];
	temp=fgets(dataline, 10000, fspecies); /* Read in the header line */
	i=0;
	while (fgets(dataline, 10000, fspecies) != NULL )
	{
		sscanf(dataline, "%s %s %d %d %le %lf %d %lf %le", (species+i)->name, (species+i)->type, &((species+i)->num), &((species+i)->mass), &((species+i)->mix), &((species+i)->upper), &((species+i)->lowertype), &((species+i)->lower), &((species+i)->lower1));
		printf("%s %s %d %d %le %lf %d %lf %le\n",(species+i)->name, (species+i)->type, (species+i)->num, (species+i)->mass, (species+i)->mix, (species+i)->upper, (species+i)->lowertype, (species+i)->lower, (species+i)->lower1);
        if (strcmp("X",species[i].type)==0) {numx=numx+1; }
		if (strcmp("F",species[i].type)==0) {numf=numf+1; }
		if (strcmp("C",species[i].type)==0) {numc=numc+1;}
		if (strcmp("A",species[i].type)==0) {numx=numx+1; numa=numa+1;}
		i=i+1;
    }
    fclose(fspecies);
	nums=numx+numf+numc;
	nn=zbin*numx;
	double Con[nn+1], fvec[nn+1], fvec1[nn+1], ConC[zbin*numc+1], Conf[zbin*numf+1]; 
	printf("%s\n", "The species list is imported.");
	printf("%s %d\n", "Number of species in model:", nums);
	printf("%s %d\n", "Number of species to be solved in full:", numx);
	printf("%s %d\n", "In which the number of aerosol species is:", numa);
	printf("%s %d\n", "Number of species to be solved in photochemical equil:", numf);
	printf("%s %d\n", "Number of species assumed to be constant:", numc);
	
	/* Variable initialization */
	int labelx[numx+1], labelc[numc+1], labelf[numf+1], MoleculeM[numx+1], listFix[numx+1], listAER[numa+1], AERCount=1; /* Standard number list of species */
	double Upflux[numx+1], Loflux[numx+1], Depo[numx+1], ConFix[numx+1], mixtemp;
	FILE *fimport;
	FILE *fimportcheck;
	
	
	/* Get the reaction list */
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
		printf("%d %s %d\n", (React+i)->dum, React[i].type, React[i].num);
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
	
	/* Load the standard reaction list and information */
	GetReaction();
	printf("%s\n", "The standard reaction databases are imported.");
	
	
	/* get the cross sections and quantum yields of molecules */   
	cross=dmatrix(1,nump,0,NLAMBDA-1);
	crosst=dmatrix(1,nump,0,NLAMBDA-1);
	qy=dmatrix(1,nump,0,NLAMBDA-1);
	qyt=dmatrix(1,nump,0,NLAMBDA-1);
    char dirroutep[1024]="Library/PhotochemOpa/";
    char photochemfile[1024];
	int stdcross[nump+1];
	double qysum[nump+1];
	fcheck=fopen("AuxillaryOut/CrossSectionCheck.dat","w");
	for (i=1; i<=nump; i++) {
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
	
	/* declare parameters before the loop */
	FILE *TPPrint;		
	double tslimit;
	
	double HenryH[NSP+1], HenryT[NSP+1], Heff;
	fhenry=fopen("Library/Henry/henry.dat","r");
	for (i=1; i<=NSP; i++) {
		fscanf(fhenry, "%lf\t", &HenryH[i]);
		fscanf(fhenry, "%lf\t", &HenryT[i]);
		printf("%e %e\n", HenryH[i], HenryT[i]);
	}
	fclose(fhenry);
	printf("%s\n","The Henry Law's constants are imported");
	
	double psH2O[zbin+1], psH2SO4[zbin+1], psS8[zbin+1];
	
	double **mixequil;
	int labels[numx+numf+1];
	
	FILE *frates;
	
	double **rad, **opt;
	
	double totalnumber, totalmass, heliumnumber;
	
	double tstep;
	double tstepold;
	int check;
	double tt, control, controlt, controlt_old, ddd, test, emax, control_nohigh;
	int controls, alindex;
	double controlz;
	double **Jaco, **Jaco1;
	int erradjust, errn, errn1, errn2, errn3;
	
	double importpres, importtemp, importmm;
	
    /* loop for different paramters */
    
    char outchemstatus[1024];
    char outoldtemp[1024];
    char crossfile[1024];
    char outstd0[1024];
    char outstd[1024];
	char dirroute[1024];
    strcpy(dirroute,OUT_DIR);
    char outphotorate[1024];
    char atomfile[1024];
    char outfile1[1024];
    char outfile2[1024];
    char outchemicalrate[1024];
    char outcolumn[1024];
    char outradiation[1024];
    char outconvergence[1024];
    char outtimescale[1024];
    char outbalance[1024];
    char outhistory[1024];
    char outnewtemp[1024];
    char outoldstd[1024];

    
    FILE *fstat;
    strcpy(outchemstatus,dirroute);
    strcat(outchemstatus,"ChemStatus.dat");
    fstat=fopen(outchemstatus,"w");
    fclose(fstat);
    
    /* Set up the P-T-z profile for calculation */
    if (TPMODE==1) {
        fp=fopen(TPLIST,"r");
        fp1=fopen(TPLIST,"r");
        s=LineNumber(fp, 1000);
        double Height[s];
        double Temp[s];
        double Pre[s];
        GetData3(fp1, 1000, s, Height, Pre, Temp);
        fclose(fp);
        fclose(fp1);
        zmin = Height[0];
        zmax = Height[s-1];
        DD=(zmax-zmin)/zbin;  /* Thickness of each layer, in unit of km */
        for (j=0; j<=zbin; j++) { z[j] = zmin+j*DD; } /* Set up the grid boundary */
        Interpolation(z, zbin+1, T, Height, Temp, s, 0);
        Interpolation(z, zbin+1, PP, Height, Pre, s, 0);
        /* make sure the both ends of the grid work well */
        T[0]=Temp[0];
        PP[0]=Pre[0];
        T[zbin]=Temp[s-1];
        PP[zbin]=Pre[s-1];
        for (i=0; i<=zbin; i++) { P[i]=pow(10,PP[i]); } /*unit: Pa*/
        /* Print-out old pressure */
        strcpy(outoldtemp,dirroute);
        strcat(outoldtemp,"/OldTemperature.dat");
        TPPrint=fopen(outoldtemp,"w");
        for (j=0; j<s; j++) {
            fprintf(TPPrint, "%f\t%f\t%f\n", Height[j], Pre[j], Temp[j]);
        }
        fclose(TPPrint);
    }
    if (TPMODE==0) {
        TPPara(P,T,TINV,zbin+1,PTOP,TTOP,PMIDDLE,TMIDDLE,PSTR,TSTR,PTROP,TTROP,PBOTTOM);
        DD=TPScale(P,T,zbin+1,z);
    }
    thickl = DD*1.0E+5; /* Thickness of each layer, in unit of cm*/
    for (j=1; j<=zbin; j++) {
        zl[j] = (z[j]+z[j-1])/2.0; /* Altitude at the center of layer */
        tl[j] = (T[j]+T[j-1])/2.0; /* Temperature at the center of layer */
        pl[j] = sqrt(P[j]*P[j-1]); /* Pressure at the center of layer */
    }
    TPPrint=fopen("AuxillaryOut/TPCheck.dat","w");
    for (j=0; j<=zbin; j++) {
        MMZ[j] = P[j]/KBOLTZMANN/T[j]*1.0E-6;
    }
    for (j=1; j<=zbin; j++) {
        MM[j]=pl[j]/KBOLTZMANN/tl[j]*1.0E-6; /*unit: Molecule cm-3*/
        printf("%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
        fprintf(TPPrint, "%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
    }
    printf("%s\n", "The Z-T-P data are imported/calculated.");
    fclose(TPPrint);
    
    /* Calculate the rainout rates throughout the atmosphere */
    for (i=1; i<=zbin; i++) {
        for (j=1; j<=NSP; j++) {
            Heff=HenryH[j]*exp(-HenryT[j]*(1/tl[i]-1/298.15)); /* temperature-dependent Henry's law constant */
            rainoutrate[i][j]=2.0E-6/55.0/NAVOGADRO/(CloudDen*1.0E-9 + 1.0/(Heff*82.05746*tl[i]))*RainF;
            /* These rates need to be multiplied by the number density of H2O to yeild a true rainout rates */
        }
    }

    /* Calculate the saturation density and condensation timescale of H2O, H2SO4 and S8 */
    waterpressure(psH2O);
    sulfuridpressure(psH2SO4);
    sulfurpressure(psS8);
    for (i=1; i<=zbin; i++) {
        nsH2O[i] = psH2O[i]/KBOLTZMANN/tl[i]*1.0E-6*SATURATIONREDUCTION;
        printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of H2O at", zl[i], "km is", psH2O[i], "Pa and", nsH2O[i], "cm-3");
    }
    for (i=1; i<=zbin; i++) {
        nsH2SO4[i] = psH2SO4[i]/KBOLTZMANN/tl[i]*1.0E-6;
        printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of H2SO4 at", zl[i], "km is", psH2SO4[i], "Pa and", nsH2SO4[i], "cm-3");
    }
    for (i=1; i<=zbin; i++) {
        nsS8[i] = psS8[i]/KBOLTZMANN/tl[i]*1.0E-6;
        printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of S8 at", zl[i], "km is", psS8[i], "Pa and", nsS8[i], "cm-3");
    }
    for (i=1; i<=zbin; i++) {
        tcondfH2SO4[i]=49.0*AMU/AERDEN*pow(8.0*KBOLTZMANN*tl[i]/PI/98.0/AMU,0.5)/AERSIZE*1.0E+6;
        tcondfS8[i]=128.0*AMU/AERDEN*pow(8.0*KBOLTZMANN*tl[i]/PI/256.0/AMU,0.5)/AERSIZE*1.0E+6;
        tcondfH2O[i]=9.0*AMU/1.0E+3*pow(8.0*KBOLTZMANN*tl[i]/PI/18.0/AMU,0.5)/AERSIZE*1.0E+6;
        printf("%s %f %s %e %s\n", "The H2SO4 condensation timescale factor at ", zl[i], "km is", tcondfH2SO4[i], "cm3 s-1");
        printf("%s %f %s %e %s\n", "The S8 condensation timescale factor at ", zl[i], "km is", tcondfS8[i], "cm3 s-1");
        printf("%s %f %s %e %s\n", "The H2O condensation timescale factor at ", zl[i], "km is", tcondfH2O[i], "cm3 s-1");
    }

    /* Import and compute Gibbs free energy of each species */
    GetGibbsForm();
    fcheckgibbs=fopen("AuxillaryOut/CheckGibbs.dat","w");
    for (i=1; i<=NSP; i++) {
        for (j=1; j<=zbin; j++) {
            fprintf(fcheckgibbs,"%2.2e\t",GibbsForm[i][j]);
        }
        fprintf(fcheckgibbs,"\n");
    }
    fclose(fcheckgibbs);
    
    /* Import the molecular opacity */
    opacCO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacSO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2O = dmatrix(1,zbin,0,NLAMBDA-1);
    opacOH = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2CO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2O2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacH2S = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacO3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCH4 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNH3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H4 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacC2H6 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacCH2O2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHCN = dmatrix(1,zbin,0,NLAMBDA-1);
    opacN2O = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNO = dmatrix(1,zbin,0,NLAMBDA-1);
    opacNO2 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacHNO3 = dmatrix(1,zbin,0,NLAMBDA-1);
    opacOCS = dmatrix(1,zbin,0,NLAMBDA-1);
    
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
    
    /* sort species-related data */
    numx1=1;
    numf1=1;
    numc1=1;
    for (i=0; i<nums; i++) {
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
            Upflux[numx1]=species[i].upper/thickl;
            Depo[numx1]=species[i].lower/thickl;
            Loflux[numx1]=species[i].lower1/thickl;
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
    
    /* Initialize composition */
    /* IMODE = 0: chemical equilibrium */
    
    if (IMODE == 0 || IMODE == 4) {
        strcpy(atomfile,dirroute);
        strcat(atomfile,"/atom.dat");
        printf("%s\t%s\n","Prepare to get atom abundance from", atomfile);
        for (i=1; i<=numx; i++) {labels[i]=labelx[i];}
        for (i=1; i<=numf; i++) {labels[numx+i]=labelf[i];}
        mixequil=dmatrix(1,zbin,1,numx+numf);
        for (j=1; j<=zbin; j++) {
            for (i=1; i<=numx+numf; i++) {
                mixequil[j][i]=0.0;
            }
        }
        chemquil(pl, tl, zbin+1, labels, numx+numf, mixequil, atomfile);
        checkmixequil(numx+numf, mixequil);
        for (j=1; j<=zbin; j++) {
            for (i=1; i<=numx; i++) {
                Con[(j-1)*numx+i]=MM[j]*mixequil[j][i];
                printf("X %d %d %e\n", j, labelx[i], mixequil[j][i]);
            }
            for (i=1; i<=numf; i++) {
                Conf[(j-1)*numf+i]=MM[j]*mixequil[j][i+numx];
                printf("F %d %d %e\n", j, labelf[i], mixequil[j][i+numx]);
            }
        }
        free_dmatrix(mixequil,1,zbin,1,numx+numf);
        Convert1(Con, ConC, Conf, labelx, labelc, labelf);
        for (i=1; i<=numx; i++) {
            if (listFix[i]==1) {
                ConFix[i]=Con[i];
            }
        }
    }
    
    /* IMODE = 4: from existing files */
    if (IMODE == 4) {
        strcpy(outstd0,dirroute);
        strcat(outstd0,"/ConcentrationSTD.dat");
        printf("%s\t%s\n","Prepare to get initial molecular concentration from", outstd0);
        fimport=fopen(outstd0, "r");
        fimportcheck=fopen("AuxillaryOut/Fimportcheck.dat","w");
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        for (j=1; j<=zbin; j++) {
            fscanf(fimport, "%lf\t", &temp1);
            fprintf(fimportcheck, "%lf\t", temp1);
            fscanf(fimport, "%lf\t", &temp1);
            fscanf(fimport, "%lf\t", &temp1);
            fscanf(fimport, "%lf\t", &importtemp);
            fscanf(fimport, "%le\t", &importpres);
            importmm = importpres/KBOLTZMANN/importtemp*1.0E-6;
            for (i=1; i<=NSP; i++) {
                fscanf(fimport, "%le\t", &xx[j][i]);
                xx[j][i]  = xx[j][i]/importmm*MM[j];
                fprintf(fimportcheck, "%e\t", xx[j][i]);
            }
            fprintf(fimportcheck,"\n");
        }
        fclose(fimport);
        fclose(fimportcheck);
        Convert2(Con, ConC, Conf, labelx, labelc, labelf);
    }
    
    strcpy(outstd,dirroute);
    strcat(outstd,"/ConcentrationSTD.dat");
    printf("%s\t%s\n","Prepare to write molecular concentration to", outstd);
    printf("%s\n", "Variable initialization successful");
    printout_std(z,outstd);
    
    strcpy(outoldstd,dirroute);
    strcat(outoldstd,"/ConcentrationSTD_Old.dat");
    printout_std(z,outoldstd);
    
    /* Mean Molecular Mass */
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
    }
    
    /* Calculate the falling velocity for aerosol species */
    fallvelocity(GA);
    for (i=1; i<=zbin; i++) {
        printf("%s %f %s %f %s\n", "The Fall velocity at", zl[i],"km is", VFall[i], "cm s-1");
    }
    
    /* Calculate the Empirical Eddy diffusion coefficient */
    if (EDDYPARA == 1) {
        for (i=1; i<zbin; i++) {
            KE[i]=KET;
            if (z[i]<ZT) {
                DenZ=MMZ[i];
            }
            if (z[i]>ZT) {
                KE[i] *= pow(DenZ/MMZ[i],0.5);
                if (KE[i]>KEH) {KE[i]=KEH;}
            }
            KE[i]=KE[i]/thickl/thickl*MMZ[i];
        }
    }
    if (EDDYPARA == 2) {
        fp=fopen(EDDYIMPORT,"r");
        fp1=fopen(EDDYIMPORT,"r");
        s=LineNumber(fp, 1000);
        double Ealt[s];
        double Etem[s];
        GetData(fp1, 1000, s, Ealt, Etem);
        fclose(fp);
        fclose(fp1);
        Interpolation(z, zbin+1, KE, Ealt, Etem, s, 0);
        for (i=1; i<zbin; i++) {
            KE[i]=KE[i]/thickl/thickl*MMZ[i];
        }
    }
    for (i=1; i<zbin; i++) {
        printf("%s %f %e\n", "Eddy diffusion at altitude", z[i], KE[i]/MMZ[i]);
    }
    /* Determine the time-step limit due to the diffusion instability */
    tslimit = thickl*thickl/KET/4.0*TSPEED;
    printf("%s %e %s\n", "The diffusion-limit stepping time is ", tslimit, "s");
    
    /* Calculate the Molecular Diffusioon coefficient */
    DM=dmatrix(1, zbin, 1, numx);
    dl=dmatrix(1, zbin, 1, numx);
    for (i=1; i<= zbin; i++) {
        for (j=1; j<=numx; j++) {
            DM[i][j]=0; /*Calculate the molecular diffusion coefficient, unit s-1*/
            if (labelx[j]==3) { /* H */
                DM[i][j]=MDIFF_H_1*1.0E+17*pow(T[i],MDIFF_H_2)/thickl/thickl;
            } /* H */
            if (labelx[j]==53) { /* H2 */
                DM[i][j]=MDIFF_H2_1*1.0E+17*pow(T[i],MDIFF_H2_2)/thickl/thickl;
            }
            dl[i][j]=DM[i][j]/2.0*((meanmolecular[i]-MoleculeM[j])*AMU*GA*thickl*1.0E-2/KBOLTZMANN/T[i]+0.38/T[i]*(tl[i+1]-tl[i]));
        }
    }
    
    /* Diffusion-limited escape velocity */
    for (j=1; j<=NSP; j++) {
        Vesc[j]=0;
    }
    for (j=1; j<=numx; j++) {
        Vesc[labelx[j]] = DM[zbin][j]/MMZ[zbin]*thickl*(meanmolecular[zbin]-MoleculeM[j])*AMU*GA*1.0E-2/KBOLTZMANN/T[zbin];
    }
    Vesc[53] = Vesc[53]*MDIFF_H2_F; /* correct for H2 photo-dissociation above the considered atmosphere */
    printf("%s %e %e %e\n", "The escape frequency of H is", Vesc[3], MMZ[zbin], T[zbin]);
    printf("%s %e\n", "The escape frequency of H2 is", Vesc[53]);
    
    /* Compute reaction rates */
    ReactionRate(); /* Calculate Reaction Rate*/
    ReactionRateM(); /* Calculate Reaction Rate*/
    ReactionRateT(); /* Calculate Reaction Rate*/
    printf("%s\n", "Forward Reaction rates initialization successful.");
    /* Compute the reverse reaction rates */
    ReverseReactionRate(zone_r,zone_m,zone_t);
    printf("%s\n", "Reverse Reaction rates initialization successful.");
    /* Print-out the reaction rates used in the model */
    frates=fopen("AuxillaryOut/ReactonRateCheck.dat", "w");
    fprintf(frates, "%s\t\t%s\t\t%s\t\t%s\n", "Type", "STD", "Rate at Bottom", "Rate at Top");
    for (i=1; i<=numr; i++) {
        fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n", "R", zone_r[i], kk[1][zone_r[i]], kk[zbin][zone_r[i]], Rkk[1][zone_r[i]], Rkk[zbin][zone_r[i]]);}
    for (i=1; i<=numm; i++) {
        fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n", "M", zone_m[i], kkM[1][zone_m[i]], kkM[zbin][zone_m[i]], RkkM[1][zone_m[i]], RkkM[zbin][zone_m[i]]);}
    for (i=1; i<=numt; i++) {
        fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n", "T", zone_t[i], kkT[1][zone_t[i]], kkT[zbin][zone_t[i]], RkkT[1][zone_t[i]], RkkT[zbin][zone_t[i]]);}
    fclose(frates);
    
    /* Generate the first radiation field */
    rad=dmatrix(0, iradmax, 0, zbin);
    opt=dmatrix(0, iradmax, 0, zbin);
    RadTransfer(rad, opt, stdcross, qysum, cross, crosst, iradmax+1);
    fout=fopen("AuxillaryOut/Radiation0.dat", "w");
    fout1=fopen("AuxillaryOut/OpticalDepth0.dat","w");
    for(ii=0; ii<=zbin; ii++)
    {
        fprintf(fout, "%s %f\n", "The initial radiation at z=", z[ii]);
        fprintf(fout1, "%s %f\n", "The optical depth at z=", z[ii]);
        for (jj=0; jj<=iradmax; jj++) {
            fprintf(fout, "%f %e\n", wavelength[jj], rad[jj][ii]); /*write the initial radiation to a file*/
            fprintf(fout1, "%f %f\n", wavelength[jj], opt[jj][ii]);} /*write the initial radiation to a file*/
    }
    fclose(fout);
    fclose(fout1);
    
    /* Compute the photolysis rates */
    JJ=dmatrix(1, zbin, 1, nump);
    Photorate(rad, cross, crosst, qy, qyt, zone_p, nump, JJ);
    printf("%s\n", "The photolysis rate is calculated for the first time.");
    strcpy(outphotorate,dirroute);
    strcat(outphotorate,"/Photorate.dat");
    printf("%s\t%s\n","Prepare to write photolysis rate to", outphotorate);
    ftemp=fopen(outphotorate, "w");
    for (ii=1; ii<=nump; ii++) {
        fprintf(ftemp, "%s\t%d\t%s\t%e\n", "The photolysis STD number", zone_p[ii], "rate at the top layer is", JJ[zbin][ii]);}
    fclose(ftemp);
    
    /* printout initial setup */
    strcpy(outfile1,dirroute);
    strcat(outfile1,"/Conx.dat");
    strcpy(outfile2,dirroute);
    strcat(outfile2,"/Conf.dat");
    fout21=fopen(outfile1,"w");
    fout22=fopen(outfile2,"w");
    fprintf(fout21,"z\n\n");
    fprintf(fout22,"z\n\n");
    fclose(fout21);
    fclose(fout22);
    strcpy(outradiation,dirroute);
    strcat(outradiation,"/Radiation.dat");
    strcpy(outcolumn,dirroute);
    strcat(outcolumn,"/ColumnDensity.dat");
    printf("%s\t%s\t%s\n","Prepare to write concentration to", outfile1, outfile2);
    printout(labelx, labelf, rad, iradmax, z, outfile1, outfile2, outradiation, outcolumn);
    strcpy(outchemicalrate,dirroute);
    strcat(outchemicalrate,"/ChemicalRate.dat");
    printoutrate(zone_r, zone_m, zone_t, JJ, zone_p, outchemicalrate);
    
    /* Main Photochemistry Loop */
    tt = 0.0;
    tstep = fmin(TSINI, tslimit);
    check = 1;
    Jaco=dmatrix(1, nn, 1, nn);
    Jaco1=dmatrix(1, nn, 1, nn);
    erradjust=0;
    if (FINE1!=1) { erradjust=numx; }
    errn=numx+1; /* choose a species that you don't want to be considered in error. set to >numx if want to consider all */
    errn1=numx+1;
    errn2=numx+1;
    errn3=numx+1;
    
    strcpy(outconvergence,dirroute);
    strcat(outconvergence,"/Convergence.dat");
    strcpy(outhistory,dirroute);
    strcat(outhistory,"/History.dat");
    fout4=fopen(outconvergence,"w");
    fout5=fopen(outhistory, "w");
    
    strcpy(outtimescale,dirroute);
    strcat(outtimescale,"/Timescale.dat");
    
    strcpy(outbalance,dirroute);
    strcat(outbalance,"/GlobalBalance.dat");
    
    for (j=0; j<NMAX; j++) {
        
        Photorate(rad, cross, crosst, qy, qyt, zone_p, nump, JJ);
        
        ChemEqu(Con, ConC, Conf, fvec1, Jaco, labelx, labelc, labelf, listAER,
                zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo, listFix);
        
        emax=1;
        tstep=2*tstep;
        
        while (emax>0.3) {
            
            /* singularflag=0;
             
             while (singularflag==0) { */
            
            tstep=0.5*tstep; /* Decrease the timestep if error larger than 0.3 */
            
            for (ii=1; ii<=nn; ii++) {
                for (jj=1; jj<=nn; jj++) {
                    Jaco1[ii][jj] = -tstep*Jaco[ii][jj];
                }
                Jaco1[ii][ii] += 1.0;
            }
            
            for (ii=1; ii<=nn; ii++) { fvec[ii]=fvec1[ii]; }
            
            singularflag=BTridiagonalm(Jaco1, fvec, zbin, numx);
            
            /*}*/
            
            /*	for (ii=1; ii<=nn; ii++) {
             printf("%d %le\n", ii, fvec[ii]);
             } */
            
            /* Calculate the time step control parameter */
            control=0;
            controlz=zmin;
            controls=1;
            controlt_old=controlt;
            controlt=0;
            for (ii=1+erradjust; ii<=nn-erradjust; ii++) {
                if ( (ii<=numx & listFix[ii]!=1) || ii>numx ) {
                    test=fabs(fvec[ii]/fmax(MINNUM, Con[ii]));
                    if (test>control & fmod(ii,numx)!=errn & fmod(ii,numx)!=errn1 & fmod(ii,numx)!=errn2 & fmod(ii,numx)!=errn3) {
                        controlt=control;
                        control=test;
                        /* if (fvec[ii]<0) {
                         controlt=test;
                         } */
                        alindex=trunc((ii-1)/numx)+1;
                        controlz=zl[alindex]; /* Get the fastest decreasing specie */
                        controls=fmod(ii, numx);
                        if (controls==0) {controls=numx;}
                    }
                    if (test<control & test>controlt & fmod(ii,numx)!=errn & fmod(ii,numx)!=errn1 & fmod(ii,numx)!=errn2 & fmod(ii,numx)!=errn3) {
                        controlt=test;
                    }
                }
            }
            if (control/controlt>100 & FINE2!=1) {
                control=controlt;
                printf("%s\n", "disregard the fastest varying point");
            }
            emax=control*tstep;
            
            /* Calculate the convergence criterium control parameter */
            control_nohigh=0;
            for (ii=1+erradjust; ii<=nn-erradjust-8*numx; ii++) {
                if ( (ii<=numx & listFix[ii]!=1) || ii>numx ) {
                    test=fabs(fvec[ii]/fmax(MINNUM, Con[ii]));
                    if (test>control_nohigh & fmod(ii,numx)!=errn & fmod(ii,numx)!=errn1 & fmod(ii,numx)!=errn2 & fmod(ii,numx)!=errn3) {
                        control_nohigh=test;
                    }
                }
            }
            
            
        }
        
        if (singularflag==0) {
            Convert1(Con, ConC, Conf, labelx, labelc, labelf);
            printf("%s\n", "encounter singular matrix");
            fprintf(fstat, "%s\n", "encounter singular matrix");
            /* record the result at each step */
            if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
                fprintf(fout5,"%f\t",tt);
                for (ii=1; ii<=zbin; ii++) {
                    for (jj=1; jj<=NSP; jj++) {
                        fprintf(fout5, "%e\t", xx[ii][jj]);
                    }
                }
                for (ii=1; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", tl[ii]);
                }
                for (ii=0; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", z[ii]);
                }
                fprintf(fout5,"\n");;
            }
            break;
        }
        
        fprintf(fout4, "%e %e %e %e %f %d\n", tt, control_nohigh, control, controlt, controlz, labelx[controls]);
        
        if (fmod(j, NPRINT)==0) { printout_timescale(Con, fvec, outtimescale);}
        
        if (j>10 & control_nohigh<Tol2) {
            check=0;
            Convert1(Con, ConC, Conf, labelx, labelc, labelf);
            printf("%s\n", "converged!");
            fprintf(fstat,"%s\n", "converged!");
            /* record the result at each step */
            if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
                fprintf(fout5,"%f\t",tt);
                for (ii=1; ii<=zbin; ii++) {
                    for (jj=1; jj<=NSP; jj++) {
                        fprintf(fout5, "%e\t", xx[ii][jj]);
                    }
                }
                for (ii=1; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", tl[ii]);
                }
                for (ii=0; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", z[ii]);
                }
                fprintf(fout5,"\n");;
            }
            break;
        }
        
        if (tt > NMAXT) {
            Convert1(Con, ConC, Conf, labelx, labelc, labelf);
            printf("%s\n", "reach the maximum integration time");
            fprintf(fstat, "%s\n", "reach the maximum integration time");
            /* record the result at each step */
            if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
                fprintf(fout5,"%f\t",tt);
                for (ii=1; ii<=zbin; ii++) {
                    for (jj=1; jj<=NSP; jj++) {
                        fprintf(fout5, "%e\t", xx[ii][jj]);
                    }
                }
                for (ii=1; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", tl[ii]);
                }
                for (ii=0; ii<=zbin; ii++) {
                    fprintf(fout5, "%f\t", z[ii]);
                }
                fprintf(fout5,"\n");;
            }
            break;
        }
        
        for (ii=1; ii<=nn; ii++) { Con[ii] = fmax(Con[ii]+tstep*fvec[ii], 0.0);} /* stepping! */
        
        /* Correct the lower boundary with fixed mixing ratio */
        /* for (ii=1; ii<=numx; ii++) {
         if (listFix[ii]==1) {
         Con[ii]=ConFix[ii];
         }
         } */
        
        if (j==0) {
            tt=tstep;
        }else {
            tt += tstep;
        }
        
        tstepold=tstep;
        if (emax>0.15) tstep=tstepold*0.9;
        if (emax>0.20) tstep=tstepold*0.7;
        if (emax<0.10) tstep=tstepold*1.1;
        if (emax<0.05) tstep=tstepold*1.3;
        if (emax<0.03) tstep=tstepold*1.5;
        if (emax<0.01) tstep=tstepold*2.0;
        if (emax<0.003) tstep=tstepold*5.0;
        if (emax<0.001) tstep=tstepold*10.0;
        
        if (tstep>TMAX) { tstep=TMAX; }
        if (tstep<TMIN) { tstep=TMIN; }
        if (tstep>tslimit) { tstep=tslimit; }
        
        Convert1(Con, ConC, Conf, labelx, labelc, labelf);
        
        /* record the result at each step */
        if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
            fprintf(fout5,"%f\t",tt);
            for (ii=1; ii<=zbin; ii++) {
                for (jj=1; jj<=NSP; jj++) {
                    fprintf(fout5, "%e\t", xx[ii][jj]);
                }
            }
            for (ii=1; ii<=zbin; ii++) {
                fprintf(fout5, "%f\t", tl[ii]);
            }
            for (ii=0; ii<=zbin; ii++) {
                fprintf(fout5, "%f\t", z[ii]);
            }
            fprintf(fout5,"\n");;
        }
        
        printf("%s %d %s %e %s %e %s\n", "finish loop", j+1, "at time", tt, "s with the timestep", tstepold, "s");
        
        RadTransfer(rad, opt, stdcross, qysum, cross, crosst, iradmax+1);
        
        if (fmod(j, NPRINT)==0) {
            printout(labelx, labelf, rad, iradmax, z,outfile1,outfile2,outradiation,outcolumn);
            printout_c();
            printout_std(z,outstd);
            printoutrate(zone_r, zone_m, zone_t, JJ, zone_p,outchemicalrate);
            printf("%s\n","print out");
            
            GlobalBalance(Con, labelx, listAER, zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo, outbalance);
            
            /* Aug 7, 2103*/
            strcpy(outnewtemp,dirroute);
            strcat(outnewtemp,"/NewTemperature.dat");
            printf("%s\t%s\n","Prepare to print out new temperature to", outnewtemp);
            GreyTemp(P,outnewtemp,TINTSET);
            
            
        }
        
    }
    
    fclose(fout4);
    fclose(fout5);
    
    /* General printout */
    printout(labelx, labelf, rad, iradmax, z,outfile1,outfile2,outradiation,outcolumn);
    printout_c();
    printout_std(z,outstd);
    /* printoutrate; */
    printoutrate(zone_r, zone_m, zone_t, JJ, zone_p,outchemicalrate);
    /* printout global balance */
    GlobalBalance(Con, labelx, listAER, zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo,outbalance);
    printf("%s\n","print out");
    
    /* Aug 7, 2013 */
    strcpy(outnewtemp,dirroute);
    strcat(outnewtemp,"/NewTemperature.dat");
    printf("%s\t%s\n","Prepare to print out new temperature to", outnewtemp);
    
    /* write out the convergence criterium */
    fprintf(fstat,"%s %e %s %e\n", "Variation is", control_nohigh, "and Duration is", tt);
    fclose(fstat);
    
    /* clean up */
    free_dmatrix(Jaco, 1, nn, 1, nn);
    free_dmatrix(Jaco1, 1, nn, 1, nn);
    free_dmatrix(JJ,1, zbin, 1,nump);
    free_dmatrix(opacCO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacSO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2O,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacOH,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2CO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2O2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacH2S,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacO3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCH4,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNH3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H4,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacC2H6,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacCH2O2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHCN,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacN2O,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNO,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacNO2,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacHNO3,1,zbin,0,NLAMBDA-1);
    free_dmatrix(opacOCS,1,zbin,0,NLAMBDA-1);
    free_dmatrix(DM, 1, zbin, 1, numx);
    free_dmatrix(dl, 1, zbin, 1, numx);
    free_dmatrix(rad, 0, iradmax, 0, zbin);
    free_dmatrix(opt, 0, iradmax, 0, zbin);
	free_dmatrix(cross,1,nump,0,NLAMBDA-1);
	free_dmatrix(qy,1,nump,0,NLAMBDA-1);
	free_dmatrix(crosst,1,nump,0,NLAMBDA-1);
	free_dmatrix(qyt,1,nump,0,NLAMBDA-1);
	
	
}
