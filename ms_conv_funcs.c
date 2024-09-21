/*******************************************************************
 * File: ms_conv_funcs.c
 * PURPOSE: compilation of main functions of convetion scheme update in 2022.
 * Original Author: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Version: v01, 2022: Original upload
 *******************************************************************
*/                  

#include <math.h>
#include "constant.h"

//=========================================================
//=== Function identifiers ================================

void ms_adiabat(int lay, double lapse[], double xxHe, double* cp); // calculate lapse rate
void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl);   // Checks each pair of layers for convective stability
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[]); // Adjust temperatures due to convection
//=== Helper function identifiers =========================
double ms_psat_h2o(double temp); //calc saturation pressure of water
double ms_psat_nh3(double temp); //calc saturation pressure
double ms_latent(int mol, double temp); //calc saturation pressure of water

//=========================================================
//=== Functions ===========================================
void ms_adiabat(int lay, double lapse[], double xxHe, double* cp)
{
// Calculate adiabat acc. to Graham+2021 Eq.(1) 
// The formulation below should be valid for dry, moist, and multi-species condensate conditions in dilute and non-dilute cases.
//
// Condensibles by species number from chemistry:
// H2O : 7
// NH3 : 9
// CO  : 20
// CH4 : 21
// CO2 : 52
// H2  : 53
// O2  : 54
// N2  : 55
// Cp avail: Air(0.8 N2, 0.2 O2), CO2, He, N2, NH3, CH4, H2, O2, CO, H2O

    //LOCAL parameters
    int i,j;
//input file    int NCONDENSIBLES = 1; //number of potential condensibles; here H2O
    double psat[NCONDENSIBLES]; //saturation pressures
    int saturated[NCONDENSIBLES]; //check which molecules are available for condensation
    double dp[NCONDENSIBLES]; //condensed fraction
    double beta[NCONDENSIBLES], latent[NCONDENSIBLES]; // L/RT for each condensible
    double cpxx; //cp * xx for all dry species for the cp_d calculation
    double MM_cp; // MM for all species we have a cp for
    double Xd,cp_d; //non-condensible gas, mol frac and heat capa.
    double Xv[NCONDENSIBLES],Xvold,cp_v[NCONDENSIBLES]; //condensibles mol fractions in vapor form and heat cap.
    double Xc[NCONDENSIBLES],cp_c[NCONDENSIBLES]; //condensibles fractions in condensed form
    double alpha[NCONDENSIBLES]; // condensate retention fraction (1 - rain-out) --- assigned with heat capacities
    double reservoirs[NCONDENSIBLES]; //keeping track of what was rained out + was available at surface at beginning

    double lapse_num, lapse_denom;
    double cp_num, cp_denom;
    double sum_beta_xv, big_sum_denom_num;
    
    //START function
    //Assign mol fractions:
    if(lay==1) printf("%s\n","--- Adiabat calculation ---");
    Xd = 1.0;
    cp_d = 0.0;
    if(lay==1) printf("%s\t%d\n","NCONDENSIBLES = ",NCONDENSIBLES);
    for (i=0; i<NCONDENSIBLES; i++)
    {
        if(lay==1) printf("%s%d%s\t%d\n","CONDENSIBLE[",i,"] = ",CONDENSIBLES[i]);

        Xc[i] = clouds[lay][CONDENSIBLES[i]]/MM[lay]; //preexisting clouds taken into account
        //Xc[i] = 0.0; //no preexisting clouds taken into account
        Xv[i] = xx[lay][CONDENSIBLES[i]]/MM[lay]; //gas fraction of condensibles

//        if(lay==1) printf("%s%d%s\t%e\n","Xv[",i,"] = ",Xv[i]);
        if(CONDENSIBLES[i]==7) psat[i] = ms_psat_h2o(tl[lay]);
        if(CONDENSIBLES[i]==9) psat[i] = ms_psat_nh3(tl[lay]);
        dp[i] = pl[lay]*Xv[i] - psat[i];
//        if(lay==1 || dp[i]>0.0) printf("%s%d%s%e\t%s%d%s%e\n","pl[",lay,"] = ",pl[lay],"psat[",i,"] = ",psat[i]);
        if(dp[i]>=0.0) 
        {
            Xc[i] += (1.0 - psat[i]/(pl[lay]*Xv[i]))*Xv[i];
            Xvold = Xv[i];
            Xv[i] *= psat[i]/(pl[lay]*Xv[i]); // here for surface layer: assumes "infinite" ocean reservoir for all condensibles
            printf("%s%d\t%s%d%s%e%s%e\t%s%d%s%e\n","Saturation in layer ",lay," Xv[",CONDENSIBLES[i],"] was ",Xvold, " now ",Xv[i]," Xc[",CONDENSIBLES[i],"] = ",Xc[i]);
        }
        else //ensure unsaturated treatment, i.e. dry adiabat, without any latent heat release is applied
        {
            Xc[i] = 0.0;
            Xv[i] = 0.0;
        }
        Xd -= Xv[i] + Xc[i];
//    if(lay==1 || dp[i]>0.0) printf("%s%e\n","Xd = ",Xd);
    }
    
    //Heat capacities, condensate retention, latent heat:
//    if(lay==1) printf("%s%e\n","xxHe = ",xxHe/MM[lay]);
    cpxx = xxHe*HeHeat(tl[lay])+xx[lay][7]*H2OHeat(tl[lay])+xx[lay][9]*NH3Heat(tl[lay])+xx[lay][20]*COHeat(tl[lay])+xx[lay][21]*CH4Heat(tl[lay])+xx[lay][52]*CO2Heat(tl[lay])+xx[lay][53]*H2Heat(tl[lay])+xx[lay][54]*O2Heat(tl[lay])+xx[lay][55]*N2Heat(tl[lay]);
    MM_cp = xxHe+xx[lay][7]+xx[lay][9]+xx[lay][20]+xx[lay][21]+xx[lay][52]+xx[lay][53]+xx[lay][54]+xx[lay][55];
    
    cp_d = cpxx / MM_cp; //this would be for dry adiabat

    for (i=0; i<NCONDENSIBLES; i++)
    {
        if(CONDENSIBLES[i]==7)  {cp_v[i] = H2OHeat(tl[lay]); cp_c[i] = H2OHeat(tl[lay]); MM_cp-=xx[lay][7];  cpxx-=xx[lay][7]*H2OHeat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==9)  {cp_v[i] = NH3Heat(tl[lay]); cp_c[i] = NH3Heat(tl[lay]); MM_cp-=xx[lay][9];  cpxx-=xx[lay][9]*NH3Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==20) {cp_v[i] = COHeat(tl[lay]);  cp_c[i] = COHeat(tl[lay]);  MM_cp-=xx[lay][20]; cpxx-=xx[lay][20]*COHeat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==21) {cp_v[i] = CH4Heat(tl[lay]); cp_c[i] = CH4Heat(tl[lay]); MM_cp-=xx[lay][21]; cpxx-=xx[lay][21]*CH4Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==52) {cp_v[i] = CO2Heat(tl[lay]); cp_c[i] = CO2Heat(tl[lay]); MM_cp-=xx[lay][52]; cpxx-=xx[lay][52]*CO2Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==53) {cp_v[i] = H2Heat(tl[lay]);  cp_c[i] = H2Heat(tl[lay]);  MM_cp-=xx[lay][53]; cpxx-=xx[lay][53]*H2Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==54) {cp_v[i] = O2Heat(tl[lay]);  cp_c[i] = O2Heat(tl[lay]);  MM_cp-=xx[lay][54]; cpxx-=xx[lay][54]*O2Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
        if(CONDENSIBLES[i]==55) {cp_v[i] = N2Heat(tl[lay]);  cp_c[i] = N2Heat(tl[lay]); MM_cp-=xx[lay][55];  cpxx-=xx[lay][55]*N2Heat(tl[lay]); 
            alpha[i] = 1.0;} //alpha assumed
    latent[i] = ms_latent(CONDENSIBLES[i],tl[lay]);
    beta[i] =  latent[i] / R_GAS / tl[lay];
    if(dp[i]>0.0) cp_d = cpxx / MM_cp; //now the non-condensing portion of moist adiabat

    //if(lay==1 || dp[i]>0.0) printf("%s%e\n","cp_d = ",cp_d);
    //if(lay==1 || dp[i]>0.0) printf("%s%e\n","cp_v = ",cp_v[0]);
    //if(lay==1 || dp[i]>0.0) printf("%s%e\n","cp_c = ",cp_c[0]);
    //if(lay==1 || dp[i]>0.0) printf("%s%d%s\t%e\n","latent[",i,"] = ",latent[i]);
    //if(lay==1 || dp[i]>0.0) printf("%s%d%s\t%e\n","alpha[",i,"] = ",alpha[i]);
    //if(lay==1 || dp[i]>0.0) printf("%s%d%s\t%e\n","beta[",i,"] = ",beta[i]);
/*    if(lay==1) printf("%s%e\n","cp_d = ",cp_d);
    if(lay==1) printf("%s%e\n","cp_v = ",cp_v[0]);
    if(lay==1) printf("%s%e\n","cp_c = ",cp_c[0]);
    if(lay==1) printf("%s%d%s\t%e\n","latent[",i,"] = ",latent[i]);
    if(lay==1) printf("%s%d%s\t%e\n","alpha[",i,"] = ",alpha[i]);
    if(lay==1) printf("%s%d%s\t%e\n","beta[",i,"] = ",beta[i]);
*/    
    }
    
    // Calculate adiabat:
    lapse_num = Xd;
    sum_beta_xv = 0.0;
    big_sum_denom_num = cp_d * Xd;
    cp_num = cp_d*Xd;
    cp_denom = Xd;

    for (i=0; i<NCONDENSIBLES; i++) 
    {
        lapse_num += Xv[i];
        sum_beta_xv += beta[i]*Xv[i];
        big_sum_denom_num += Xv[i]*(cp_v[i] - R_GAS*beta[i] + R_GAS*beta[i]*beta[i]) + alpha[i]*Xc[i]*cp_c[i];
        cp_num +=  Xv[i]*cp_v[i] + alpha[i]*Xc[i]*cp_c[i];
        cp_denom += Xv[i];
    }

    lapse_denom = Xd * big_sum_denom_num / (R_GAS*(Xd + sum_beta_xv)) + sum_beta_xv;

    lapse[lay] = lapse_num / lapse_denom;

    *cp = cp_num / cp_denom; //return to climate module

    //if(lay==1) printf("%s%e\n","lapse_num = ",lapse_num);
    //if(lay==1) printf("%s%e\n","lapse_denom = ",lapse_denom);
    //if(lay==1) printf("%s%e\n","sum_beta_xv = ",sum_beta_xv);
    //if(lay==1) printf("%s%e\n","big_sum_denom_num = ",big_sum_denom_num);
//    if(lay==1) printf("%s%e\n","lapse = ",lapse[lay]);

    //propagate compositional changes back to main program 
    for (i=0; i<NCONDENSIBLES; i++)
    {
        if(dp[i]>0.0)
        {
            xx[lay][CONDENSIBLES[i]] = Xv[i] * MM[lay]; //since xv set to 0.0 if unsaturated
        }
        clouds[lay][CONDENSIBLES[i]] = alpha[i] * Xc[i] * MM[lay]; //if unsaturated, then clouds should dissolve?
    }

    //Now correct all xx for any rained out molecules:
    //removal of condensate through rain out would cause SUM(xx) /= 1.0 otherwise
    //(or to date 07/2022 it would actually cause heliumfraction in ms_rad_conv to increase by (1-alpha)*Xc)
    for (j=0; j<NSP+1; j++)
    {
        for (i=0; i<NCONDENSIBLES; i++)
        {
            xx[lay][j] = xx[lay][j] / (1 - ( 1-alpha[i])*Xc[i]); 
        }
        // more consistent would be to adjust MM down - and Presuure with it. But careful with scenarios where Psurf is a function of psat*RH dependent ocean evaporation
    }

//    printf("%s %d %s %d\n","ms_ncl",*ncl,"ms_nrl",*nrl);
}// END: void ms_adiabat()
//*********************************************************
//*********************************************************
void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl)
{
// Check each pair of layers for convective stability
//ncl, nrl, and isconv are handed over as 0 initially

    //LOCAL parameters
    int i;
    for (i=0; i<=zbin; i++) isconv[i] = 0;
    
    //TESTING!!
    //isconv[0] = 1;

    //START function
    for (i=0; i<zbin; i++)
    {
//        printf("%s %d%s %f%s%f\n","Layer",i,": Tabove=",tempb[i+1]," but adiabatically ",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
//        if( tempb[i+1] < ( tempb[i] * pow(P[i+1]/P[i],lapse[i+1]) /0.98 ) )
        if( tempb[i+1] < ( tempb[i] * pow(P[i+1]/P[i],lapse[i+1]) +1e-1) ) //for numerical stability
        {
        //printf("%s %d\t%s %d\t%f %s %f\n","Layer",i+1,"isconv",isconv[i],tempb[i+1],"<",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
            isconv[i+1] = 1;
            //isconv[i] = 1;
        }
//    printf("%s %d %s %d\n","ms_i",i,"ms_isconv",isconv[i]);
    }
    for (i=1;i<=zbin;i++)
    {
        if(isconv[i]) *ncl+=1;
        else *nrl+=1;
    }
//    printf("%s %d %s %d\n","ms_ncl",*ncl,"ms_nrl",*nrl);
}// END: void ms_conv_check()
//*********************************************************
//*********************************************************
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[])
{
// Adjust temperatures due to convection

    //LOCAL parameters
//    int ncreg[zbin+1] = {0}; //determine convective regions
    int i,j,cr,nccount = 0;
    double potT[zbin+1]={0};//, pot_temp[zbin+1];
    double enthalpy[zbin+1]={0}, potTint[zbin+1]={0};
    double ppk[zbin+1];
    double t_old[zbin+1];

    //START function
    for (i=0;i<=zbin;i++)
    {   
        t_old[i]=tempb[i]; //copy for later
        ppk[i]=1.0;
    }

//    printf("%s\t%s\n","i","ppk[i]");
    if(isconv[0]) {nccount+=1; ncreg[0] = nccount;}
    for (i=1; i<=zbin; i++) //number convective regions
    {
//        printf("%d\t%e\n",i,ppk[i]);
        if(isconv[i]==1) 
        {
            if(isconv[i-1]==0) nccount += 1;
            ncreg[i] = nccount;
        }
    }//END: number convective regions
//    print:f("%s\t%s\t%s\n","Layer","isconv","Conv_Region");
//    for (i=zbin; i>=1; i--) printf("%d\t%d\t%d\n",i,isconv[i],ncreg[i]);

    //CALC potential Temps
//    printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

    for (i=1; i<=zbin; i++)
    {
//temps too high:        enthalpy[ncreg[i]] += cp[i] * (tempb[i]+tempb[i-1])/2.0 * (P[i-1] - P[i]);
        enthalpy[ncreg[i]] += cp[i] * tempb[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-2","potT enthalpy integral-2");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

        j=i;
//        printf("%s\t%s\t%s\t%s\t%s\t%s\n","i","j","ppk","Pj","Pj-1","lapse");
        while (abs(ncreg[j]-ncreg[i]) <1.0 && isconv[i]==1 && j>0)
        {
            ppk[i] *= pow(P[j] / P[j-1],lapse[j]);
//            printf("%d\t%d\t%e\t%e\t%e\t%f\n",i,j,ppk[i],P[j],P[j-1],lapse[j]);
            j -= 1;
        }
        potTint[ncreg[i]] += cp[i] * ppk[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-3","potT enthalpy integral-3");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}
        
//        printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//        for (cr=0; cr<=nccount; cr++) printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);
    }
    
//    printf("%s\t%s\t%s\t%s\n","region","enthalpy","potT enthalpy integral","potT");

    for (cr=1; cr<=nccount; cr++)
    {
        potT[cr] = enthalpy[cr] / potTint[cr];
//        printf("%d\t%e\t%e\t%e\n",cr,enthalpy[cr],potTint[cr],potT[cr]);
    }

    //ADJUST temps
    for (i=1; i<=zbin; i++) //first mid-layer:
    {
        pot_temp[i] = potT[ncreg[i]];
        if(isconv[i]) 
        {
            tempb[i] = potT[ncreg[i]] * ppk[i];
            if (isconv[i-1]==0) tempb[i-1] = potT[ncreg[i]];
        }
    }
/*    for (i=1; i<zbin; i++) //now layer boundaries
    {
        tempb[i] = 0.5 * ( tl[i] + tl[i+1]);
    }*/
    //TOA:
    //tempb[zbin] = tl[zbin] + (tl[zbin] - tempb[zbin-1]);
    //ADJUST surface
    //if(isconv[1]) tempb[0] = potT[ncreg[1]] * (1 + 0.5 * (pow(P[0]/(P[1]),lapse[1]) - 1.0));
    // covered already above if(isconv[1]) tempb[0] = potT[ncreg[1]];

    printf("%s\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\n","Layer","cp","lapse","isconv","Conv_Region","Pot_Temp","T_old","T_new");
    for (i=zbin; i>=0; i--) {printf("%d\t%f\t%f\t%d\t%d\t\t%e\t%f\t%f\n",i,cp[i],lapse[i],isconv[i],ncreg[i], pot_temp[i], t_old[i], tempb[i]);}

//atexit(pexit);exit(0); //ms debugging mode
}//END: void ms_temp_adj()
//*********************************************************
//*********************************************************

//=========================================================
//=== Helper Functions ====================================
double ms_psat_h2o(double temp)
{
    double a,P;
    
    //printf("%s%f\n","tl = ",tl);
    if (temp<273.16) 
    {
	// over ice
	// Formulation from Murphy & Koop (2005)
	P = exp(9.550426-5723.265/temp+3.53068*log(temp)-0.00728332*temp); // [Pa] 
    }
    //else if(tl<373.15) 
    else
    {
	// over water
	// Formulation from Seinfeld & Pandis (2006)
	a = 1-373.15/temp;
	P = 101325*exp(13.3185*a-1.97*a*a-0.6445*a*a*a-0.1229*a*a*a*a); // [Pa]
    }
/*    else 
    {
        P = pow(10,8.14019-1810.94/(244.485+tl-273.15))*133.322387415;
    }
*/	
    //printf("%s%f\n","Psat = ",*P);
    return P;
}// END: double ms_psat_h2o()
//*********************************************************
//*********************************************************
double ms_psat_nh3(double temp)
{
    double P;
    
    //printf("%s%f\n","tl = ",tl);
    if (temp<195.40) 
    {
	// solid
	// Formulation from Lodders, Fegley (1998)
	P = pow(10.0,6.9-1588/temp) * 1.e+5; // [Pa] 
    }
    else if(temp<300.00) 
    {
        // liquid
	P = pow(10.0,5.201-1248/temp) * 1.e+5; // [Pa] 
    }
    else
    {
	// undefined - set for unsaturated case
	P = 1e+20; // [Pa]
    }
    //printf("%s%f\n","Psat = ",*P);
    return P;
}// END: double ms_psat_nh3()
//*********************************************************
//*********************************************************
double ms_latent(int mol, double temp)
{
    double l;
    double t_crit,t_triple;  //Temps [K] of boiling point, triple point, and critical point (not used so far)
    double lvap_bp, lvap_tp; //latent heats [J/kg] for liquid-vapor phase change at boiling point and triple point
    double lfuse;            //latent heat [J/kg] for liquid-solid phase change
    double lsubl;            // latent heat [J/kg] for solid-vapor phase change
    double molweight;        // g/mol of each molecule to convert to J/mol
    
// Values from Pierrehumbert 2010 2.3.2, Table2.1, unless stated otherwise:
// NIST lvap and lsubl calculated from dH[kj/mol] / MolWeight [g/mol] *1e6
//H2O
    if(mol==7)  {t_crit= 647.1; t_triple=273.16; 
        lvap_bp=22.55e+05; lvap_tp=24.93e+05; 
        lfuse=3.34e+05; lsubl=28.4e+05;
        molweight = 18.02;
    }
//NH3
    if(mol==9)  {t_crit= 405.5; t_triple=195.4; 
        lvap_bp=13.71e+05; lvap_tp=16.58e+05; 
        lfuse=3.314e+05; lsubl=19.89e+05;
        molweight = 17.03;
    } 
//CO --- Temps: NIST Chemistry Webbook SRD 69
    if(mol==20) {t_crit= 134.45; t_triple=68.12; 
        lvap_bp=2.14e+05; lvap_tp=2.14e+05 /*NIST dH=6.0 for both*/; 
        lfuse=2.98e+04 /*engineeringtoolbox.com - calculated from BTU/lb*/; 
        lsubl=2.89e+05 /*NIST dH=8.1*/;
        molweight = 28.01;
            //since lsubl=lfuse+lvap_bp+cp_liquid*(Tboil=81.65 - Tmelt=68.12) -> cp_CO_liquid =~ 4006 [J/kg/K] = 0.1122 [kJ/mol/K] ?
    } 
//CH4
    if(mol==21) {t_crit= 190.44; t_triple=90.67; 
        lvap_bp=5.1e+05; lvap_tp=5.36e+05; 
        lfuse=5.868e+04; lsubl=5.95e+05;
        molweight = 16.04;
    } 
//CO2
    if(mol==52) {t_crit= 304.2; t_triple=216.54; 
        lvap_bp=3.79e+05 /*NIST from dH=16.7*/; lvap_tp=3.97e+05; 
        lfuse=1.96e+05; lsubl=5.93e+05;
        molweight = 44.01;
    } 
//H2
    if(mol==53) {t_crit= 33.2; t_triple=13.95; 
        lvap_bp=4.54e+05; lvap_tp=4.54e+05 /*from bp since no values avail.*/; 
        lfuse=5.82e+04; 
        lsubl=5.65e+05; //=lfuse+lvap_bp +dH_liquid [=cp_liquid*(Tboil=20.39K - Ttriple)] - since no data avail.
        molweight = 2.02;
            //cp[tp]=~6.57e3[J/kg/K], cp[bp]=~9.77e3[J/Kg/K] -> dH_liquid=0.5*(6.57e3+9.77e3)*(20.39-13.95) = 5.26e4
            //lsubl=lfuse+dH_liquid+lvap_bp
            //Liquid H2 data from KAERI/TR-2723/2004 (https://inis.iaea.org/collection/NCLCollectionStore/_Public/36/045/36045728.pdf)
    } 
//O2
    if(mol==54) {t_crit= 154.54; t_triple=54.3; 
        lvap_bp=2.13e+05; lvap_tp=2.42e+05; 
        lfuse=13.9e+04; lsubl=2.56e+05;
        molweight = 32.00;
    } 
//N2
    if(mol==55) {t_crit= 126.2; t_triple=63.14; 
        lvap_bp=1.98e+05; lvap_tp=2.18e+05; 
        lfuse=2.573e+04; lsubl=2.437e+05;
        molweight = 28.01;
    } 

//Now calculate the appropriate l:
    if (temp<=t_triple) {l = lsubl;}
    else if (temp<t_crit) {l = (t_crit-temp)/(t_crit-t_triple)*lvap_tp + (temp-t_triple)/(t_crit-t_triple)*lvap_bp;}
    else if (temp>=t_crit) {l = lvap_bp;}

    l = l * molweight / 1.0e+3; // convert J/kg to J/mol
    //printf("%s%f\n","Psat = ",*P);
    return l;
}// END: couble ms_latent()
//*********************************************************
//*********************************************************

//atexit(pexit);exit(0); //ms debugging mode
//printf("%s\n", "HERE"); //template
