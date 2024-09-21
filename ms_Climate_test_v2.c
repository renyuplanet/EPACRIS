/*Function to calculate atmnospheric temperatures and condensate states in Radiative-Convective Equilibrium
 *
 * Original Author: Renyu Hu (renyu.hu@jpl.nasa.gov)
 * Version: 2021, 2-Stream update to Heng+2019 flux equations including direct stellar beam
 * Editor v2021: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Version: 2022, convection update including condensibles and new adiabat formulation from Graham+2021 (Pierrehumbert's group)
 * Editor v2022: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * **********************************************************************************************************************************
 * Current Version: 2023, 
 *                  complete restructuring of Rad-Conv loop towards equilibrium; 
 *                  adding numerous switches for different studies;
 *                  self-consistent ocean-atmosphere boundary based on Psat & RH
 * Editor v2023: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * **********************************************************************************************************************************
*/

#include <math.h>
#include "constant.h"
#include "ludcmp.c"
#include "lubksb.c"
#include "BTridiagonal.c"
#include "ms_radtrans_test.c"
#include "ms_conv_funcs.c"

void ms_Climate(double tempeq[], double P[], double T[], double Tint, char outnewtemp[],char outrtdiag[],char outrcdiag[],char outcondiag[]);

void ms_Climate(double tempeq[], double P[], double T[], double Tint, char outnewtemp[],char outrtdiag[],char outrcdiag[],char outcondiag[])
{
    int i,j,k,radcount,jabove,iconden;
    
    /* Initial profile */
    double tempb[zbin+1], tempbnew[zbin+1], tempc[zbin+1];
    for (j=0; j<=zbin; j++) {
		tempb[j] = T[j]; //ms23: double grid
		//tempb[j] = Tdoub[2*j]; //bottom-up!
	}

   /* Helium calculation */ 
    double heliumnumber[zbin+1];
    for (j=1; j<=zbin; j++) {
        heliumnumber[j] = MM[j];
        for (i=1; i<=NSP; i++) {
            heliumnumber[j] -= xx[j][i];
        }
        if (heliumnumber[j]<0.0) {
            heliumnumber[j] = 0.0;
        }
    }
    
    /* Define variables for iteration */
    double znew[zbin+1], scaleheight, GA;
    GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET;
    double lapse[zbin+1], cp[zbin+1];
    //double Rflux[zbin+1];
    double *Rflux;
    double Rfluxmax,dRfluxmax;
    int isconv[zbin+1], ncl=0, nrl=zbin; //convective layers setup
    int isconv_old[zbin+1], prevconv[zbin+1];
    int isconv_sum, deltaconv;
    double t_old[zbin+1],tl_old[zbin+1]; //for convective adj test
    int isequil[zbin+1], sumisequil;
    double dt[zbin+2]; //time step if TIME_STEPPING
    char RTstopfile[1024]; 
    int ncreg[zbin+1] = {0}; //convective regions
    double pot_temp[zbin+1] = {0.0}; //potential temperature per convective region
    int nclouds[zbin+1] = {0}; //number of condensed species per layer
    strcpy(RTstopfile,OUT_DIR);
    strcat(RTstopfile,"rt.stop");
    printf("%s %s\n", "RTstopfile = ", RTstopfile);
//atexit(pexit);exit(0); //ms debugging mode

    for (j=1;j<=zbin+1;j++) dt[j] = 1.0e+0; //for TIME_STEPPING: init in some random unit to ...
    for (j=1; j<=zbin; j++) isconv[j] = 0;
//    for (j=1; j<=zbin; j++) ms_adiabat(j,lapse,heliumnumber[j],&cp[j]); //we only need initial cp[] here if dt depends on it

    Rflux = dvector(0,nrl);

//*************************************************************
// main iterative Rad-Conv loop:
//*************************************************************
    for (i=1; i<=NMAX_RC; i++) {
        
        printf("%s\n",filleq);
        printf("%s %d\n", "Rad-Conv loop i = ", i);
        printf("%s\n",filleq);
        
        /* store last RC boudary */
        for (j=0; j<=zbin; j++) t_old[j] = tempb[j];
        for (j=1; j<=zbin; j++) {
            //tl_old[j] = tl[j]; //ms23: seems weird to overwrite here
            isconv_old[j] = isconv[j];
            prevconv[j] = isconv[j]; //for convective adjustment iteration
            isconv[j] = 0; //Allways run RT on full grid
        }
        ncl=0;nrl=zbin; //Allways run RT in full grid

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Radiative Transfer iteration
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        int pevery=20, pcount=0;//Added for printing every few steps
        double Tsaved[zbin+1];
        int saveevery=20; //save Told every n steps for dt progression
        int RTsteplimit;

        for (j=0;j<=nrl;j++) Rflux[j] = 0.0; 
        
        RTstepcount=0; //initialized in main
        if (i==1) RTsteplimit = NMAX_RT;
        if (i>1) RTsteplimit = NRT_RC;

        for (j=0;j<=zbin;j++) Tsaved[j] = tempb[j]; //store T for dt progression
        for (j=0;j<=zbin;j++) isequil[j] = 0; //store T for dt progression

        //while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || RTstepcount==0 || sumisequil < zbin) && RTstepcount<RTsteplimit && access(RTstopfile, F_OK) !=0 ) { 
        //sumisequil is not defined yet for jacobian; also not a reliable proxi yet for convergeance in time stepping
        while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || RTstepcount==0) && RTstepcount<RTsteplimit && access(RTstopfile, F_OK) !=0 ) { 
        //while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || RTstepcount==0 ) && rtstepcount<1000 ) { 
        //while (rcount<=NMAX_RC ) { 
            RTstepcount +=1;
            pcount +=1;
            sumisequil=0;

            printf("%s %d\n", "Radiative loop RTstepCount = ", RTstepcount);
            ms_RadTrans(Rflux, tempbnew, P, ncl, isconv, lapse, tempb, Tint, cp, dt, isequil);
            printf("%s\t%s\t%s\t%s\t%s\t%s\n","==== i","P[:]","T_old[:]","T_new[:]","Layer_equil","dt[:] ====");
            for (j=zbin; j>=0; j--) {
                printf("%d\t%e\t%f\t%f\t%d\t%e\n",j,P[j],tempb[j],tempbnew[j],isequil[j],dt[zbin+1-j]);
                sumisequil += isequil[j];
            }
            
            if (RTstepcount%saveevery==1) {//store T every n steps for dt progression
                for (j=0;j<=zbin;j++) {
                    //if(isequil[j]==0){//only advance dt for non-equilibrium layers
                        if (fabs(tempbnew[j] - Tsaved[j]) < 0.5*saveevery*(fabs(tempbnew[j]-tempb[j]))) {
                            dt[zbin+1-j] /=1.5; //as in HELIOS
                        } else {
                            dt[zbin+1-j] *= 1.1; //as in HELIOS
                        }
                    //}
                    Tsaved[j] = tempb[j];
                }
            }
            for (j=0; j<=zbin; j++) tempb[j]=tempbnew[j];
            for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);

            Rfluxmax=0.0;
            dRfluxmax=0.0;
            radcount = 0;
            for (k=1; k<=zbin; k++) {
                if (isconv[k] == 0) {
                    radcount += 1;
                    Rfluxmax = fmax(Rfluxmax,fabs(Rflux[radcount]));
                    dRfluxmax = fmax(dRfluxmax,fabs((Rflux[radcount-1] - Rflux[radcount])));
                }
            }
            printf("%s\t%s\t%s\t%s\n","Conv_loop","Rad_loop","Residual flux max","dFnet max");
            printf("%d\t\t%d\t\t%f\t%f\n", i, RTstepcount, Rfluxmax, dRfluxmax);
            printf("%s\n",filleq);
            if ( (Rfluxmax < Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax < Tol_RC) || (dRfluxmax < 0.2*Tol_RC)  ) {
                //if(sumisequil>=zbin){ //ms23: needs fine tuning before it can be used
                    printf("Radiative transfer done\n");
                    break;
               // }
            }

            if (pcount = pevery)
            {
                /* Print-out new TP profile*/
                //ms23: also print diagnostics
                FILE *fp,*frt;
                fp=fopen(outnewtemp,"w");
                frt=fopen(outrtdiag,"w");
                fprintf(frt,"%s\t%8s\t%8s\t%10s\t%10s\t%8s\t%s\n","Layer","Height","log(P)","Temp_bound.","Rflux","isequil","dt");
                for (j=0; j<=zbin; j++) 
                {
                    fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempb[j]); //temperature file
                    fprintf(frt, "%d\t%8f\t%8f\t%10f\t%10f\t%d\t%e\n",j, znew[j], log10(P[j]), tempb[j], Rflux[j], isequil[j], dt[j]); //RadTrans diagnostics
                }
                fclose(fp);
                fclose(frt);
                //pcount = 0; nneded for convection later
            }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } //ms22: END of radiative transfer iteration
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if(access(RTstopfile, F_OK )==0) printf("\n\n%s\n\n","===== RTstopfile found. RT loop aborted! =====");

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //ms22: convective adjustment iteration:
    //needed to avoid temperature jumps at upper/lower tails of convective regions 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        deltaconv = zbin; //for convective adjustment iteration
        while (deltaconv > 0) 
        {
            for (j=1; j<=zbin; j++) {pot_temp[j] = 0.0; ncreg[j] = 0;} 
            /* heat capacity and laspe rate */
            for (j=1; j<=zbin; j++) 
            {
                ms_adiabat(j,lapse,heliumnumber[j],&cp[j]); //calculate appropriate lapse rate depending on dry, moist, and condensate conditions
                
                //need to propagate cold traps upwards
                if (j<zbin) {
                    for (jabove=j+1; jabove<=zbin; jabove++) {
                        for (iconden=0; iconden<NCONDENSIBLES; iconden++) {
                            if (xx[jabove][CONDENSIBLES[iconden]]/MM[jabove] > xx[j][CONDENSIBLES[iconden]]/MM[j]) {
                                xx[jabove][CONDENSIBLES[iconden]] = xx[j][CONDENSIBLES[iconden]]/MM[j]*/MM[jabove];
                            }
                        }
                    }
                }
                
                //check for clouds
                nclouds[j] = 0;
                for (k=1;k<=NSP;k++)
                {
                    if(clouds[j][k] > 0.0) nclouds[j] += 1;
                }
            }

            /* new RC boundary */
            ncl = 0;
            nrl = 0;
            /* determine convection, and record convective layers */

            if(CONVEC_ADJUST==0)
            {
                for (j=zbin; j>=1; j--) {
                    //if ( tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.999999) ) {
                    if ( i!=1 && (tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.98)) ) {
                    //if ( tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.95) ) {
                        tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]); //ms2022: moved down to include manually modded layers
                        isconv[j] = 1;
                        ncl = ncl+1;
                    } else {
                        isconv[j] = 0;
                        nrl = nrl+1;
                    }
                }
            } //END: CONVEC_ADJUST=0

            if(CONVEC_ADJUST==1) ms_conv_check(tempb, P, lapse, isconv, &ncl, &nrl); //ms22: check convective layers

//          printf("%s %d %s %d\n", "ncl", ncl, "nrl", nrl);

        // remove single radiative layer between convective layers
/*          for (j=2; j<zbin; j++) { //rh->ms 2021
                if (isconv[j]==0 && (isconv[j+1]==1 && isconv[j-1]==1)) {
                    isconv[j]=1;
                    ncl = ncl+1;
                    nrl = nrl-1;
                }
            }
*/
        
        // remove single convective layer between radiative layers
/*          for (j=2; j<zbin; j++) { //rh->ms 2021
                if (isconv[j]==1 && (isconv[j+1]==0 && isconv[j-1]==0)) {
                    isconv[j]=0;
                    ncl = ncl-1;
                    nrl = nrl+1;
                }
            }
*/        
   
            if(CONVEC_ADJUST==0)
            {
                printf("%s\t%s\t\t%s\t\t%s\t%s\t%s\n","Layer","cp","lapse","isconv","T_old","T_new");
                for (j=zbin; j>=0; j--) {printf("%d\t%f\t%f\t%d\t%f\t%f\n",j,cp[j],lapse[j],isconv[j], t_old[j], tempb[j]);}
            }

            if(CONVEC_ADJUST==1 && ncl>=0) ms_temp_adj(tempb, P, lapse, isconv, cp, ncreg, pot_temp); //ms22: assign convective regimes, adjust temperatures

            //ms22: check if convective regions were fully addressed:
            deltaconv = 0;
            for (j=1;j<=zbin;j++)
            {
                if (abs(prevconv[j]-isconv[j])>0) deltaconv += 1;
                prevconv[j] = isconv[j];
            }
            for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);
            printf("%s %d\n","deltaconv=",deltaconv);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } //ms22: END of convective adjustment iteration
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        /* hydrostatic equilibrium */
        znew[0]=0.0;
        for (j=1; j<=zbin; j++) {
            tempc[j] = (tempb[j]+tempb[j-1])/2.0;
            scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000 ;
            znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);
            //ms23: double grid
            Tdoub[2*j] = tempb[j];
            Tdoub[2*j-1] = tempc[j];
        }
        Tdoub[0] = tempb[0];
        printf("%s %d %s %d\n", "conv layer", ncl, "rad layer", nrl);
        /* compare old and new RC boundary */
        isconv_sum=0;
        for (j=1; j<=zbin; j++) {
            isconv_sum += abs(isconv_old[j]-isconv[j]);
        }
        printf("%s %d\n", "change in convective layer is", isconv_sum);

        if (pcount = pevery)
        {
        //ms23: print convection diagnostics
            FILE *frc,*fcon;
            frc=fopen(outrcdiag,"w");
            fprintf(frc,"%s\t%10s\t%8s\t%8s\t%5s\t%s\t%s\t%s\t%10s\t%s\n","Layer","Height","log(P)","cp","lapse","nClouds","isconv","Conv_Region","Pot_Temp","Temp_center");
            fcon=fopen(outcondiag,"w");
            fprintf(fcon,"%s\t%8s\t%s","Layer","log(P)","nClouds");
            for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%s %3d\t%s %4d","Molec.",CONDENSIBLES[k],"Cloud",CONDENSIBLES[k]); //Condensibles diagnostics
            fprintf(fcon,"\n");

            for (j=1; j<=zbin; j++) 
            {
                fprintf(frc,"%d\t%10f\t%8f\t%8f\t%8f\t%d\t%d\t%d\t%10e\t%f\n",j,znew[j],log10(P[j]),cp[j],lapse[j],nclouds[j],isconv[j],ncreg[j],pot_temp[j],tempc[j]); //Rad-Conv diagnostics
                fprintf(fcon,"%d\t%f\t%d",j,log10(P[j]),nclouds[j]); //Condensibles diagnostics
                for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%.4e\t%.4e",xx[j][CONDENSIBLES[k]]/MM[j],clouds[j][CONDENSIBLES[k]]/MM[j]); //Condensibles diagnostics
                fprintf(fcon,"\n");
            }
            fclose(frc);
            fclose(fcon);
            pcount = 0;
        }

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //ms22: smoothing temperature profile
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for (j=1; j<zbin; j++) tempeq[j] = tempb[j-1]/3. + tempb[j]/3. + tempb[j+1]/3.;
        //ms23: different smoothing for TOA and BOA
        tempeq[0] = tempeq[1] + 1.0/3.0 * (tempeq[1] - tempeq[4]); //trend of smoothed T 1 through T 4
        tempeq[zbin] = tempeq[zbin-1] + 1.0/3.0 * (tempeq[zbin-1] - tempeq[zbin-4]); //trend of smoothed T N-1 through T N-4
        //tempeq[0] = tempb[0]/2. + tempb[1]/2.; //surf
        //tempeq[zbin] = tempb[zbin-1]/2. + tempb[zbin]/2.; //TOA
        
        for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // DETERMINE IF CONVERGED
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (isconv_sum == 0 && i!=1 && sumisequil>=zbin) {
            printf("%s\n",filleq);
            printf("Climate converged\n");
            printf("%s\n",filleq);
            break;
        }
    
        if (ncl == 0) {
            printf("%s\n",filleq);
            printf("No convective layers found. Rad-Conv-loop done, but not all layers meet selected radiative equilibrium requirement\n");
            printf("%s\n",filleq);
            break;
        }
//*************************************************************
    }// END main iterative Rad-Conv loop:
//*************************************************************
    
    /* Print-out new TP profile */
    //ms23: also diagnostic files
    FILE *fp,*frt,*frc,*fcon;
    fp=fopen(outnewtemp,"w");
    frt=fopen(outrtdiag,"w");
    fprintf(frt,"%s\t%8s\t%8s\t%10s\t%10s\t%8s\t%s\n","Layer","Height","log(P)","Temp_bound.","Rflux","isequil","dt");
    frc=fopen(outrcdiag,"w");
    fprintf(frc,"%s\t%10s\t%8s\t%8s\t%8s\t%s\t%s\t%s\t%10s\t%s\n","Layer","Height","log(P)","cp","lapse","nClouds","isconv","Conv_Region","Pot_Temp","Temp_center");
    fcon=fopen(outcondiag,"w");
    fprintf(fcon,"%s\t%8s\t%s","Layer","log(P)","nClouds");
    for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%s %3d\t%s %4d","Molec.",CONDENSIBLES[k],"Cloud",CONDENSIBLES[k]); //Condensibles diagnostics
    fprintf(fcon,"\n");

    for (j=0; j<=zbin; j++) 
    {
        fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempeq[j]);
        fprintf(frt, "%d\t%8f\t%8f\t%10f\t%10f\t%d\t%e\n",j, znew[j], log10(P[j]), tempb[j], Rflux[j], isequil[j], dt[j]); //RadTrans diagnostics
        if (j>0)
        {
            fprintf(frc,"%d\t%10f\t%8f\t%8f\t%5f\t%d\t%d\t%d\t%10e\t%f\n",j,znew[j],log10(P[j]),cp[j],lapse[j],nclouds[j],isconv[j],ncreg[j],pot_temp[j],tempc[j]); //Rad-Conv diagnostics
            fprintf(fcon,"%d\t%f\t%d",j,log10(P[j]),nclouds[j]); //Condensibles diagnostics
            for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%.4e\t%.4e",xx[j][CONDENSIBLES[k]]/MM[j],clouds[j][CONDENSIBLES[k]]/MM[j]); //Condensibles diagnostics
            fprintf(fcon,"\n");
        }
    }
    fclose(fp);
    fclose(frt);
    fclose(frc);
    fclose(fcon);
    
    free_dvector(Rflux,0,nrl);
}
//atexit(pexit);exit(0); //ms debugging mode
