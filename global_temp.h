/*----------------------- global.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: July 20, 2011

--------------------------------------------------------------------- */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

/*---- External Variables ------------------------------------------- */

extern double meanmolecular[];
extern double zl[];
extern double pl[];
extern double tl[];
extern double MM[], MMZ[];
extern double wavelength[];
extern double solar[];
extern double crossr[], crossa[3][NLAMBDA], sinab[3][NLAMBDA], asym[3][NLAMBDA];
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO; 
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4; 
extern double **opacNH3;
extern double **opacC2H2, **opacC2H4, **opacC2H6, **opacHCN, **opacCH2O2, **opacHNO3;
extern double **opacN2O, **opacN2, **opacNO, **opacNO2, **opacOCS;
extern double **opacHF, **opacHCl, **opacHBr, **opacHI, **opacClO, **opacHClO;
extern double **opacHBrO, **opacPH3, **opacCH3Cl, **opacCH3Br, **opacDMS, **opacCS2;
extern double MeanCO2[zbin+1], MeanO2[zbin+1], MeanSO2[zbin+1], MeanH2O[zbin+1], MeanOH[zbin+1], MeanH2CO[zbin+1];
extern double MeanH2O2[zbin+1], MeanHO2[zbin+1], MeanH2S[zbin+1], MeanCO[zbin+1], MeanO3[zbin+1], MeanCH4[zbin+1];
extern double MeanNH3[zbin+1];	
extern double MeanC2H2[zbin+1], MeanC2H4[zbin+1], MeanC2H6[zbin+1], MeanHCN[zbin+1], MeanCH2O2[zbin+1], MeanHNO3[zbin+1];
extern double MeanN2O[zbin+1], MeanN2[zbin+1], MeanNO[zbin+1], MeanNO2[zbin+1], MeanOCS[zbin+1];
extern double SMeanCO2[zbin+1], SMeanO2[zbin+1], SMeanSO2[zbin+1], SMeanH2O[zbin+1], SMeanOH[zbin+1], SMeanH2CO[zbin+1];
extern double SMeanH2O2[zbin+1], SMeanHO2[zbin+1], SMeanH2S[zbin+1], SMeanCO[zbin+1], SMeanO3[zbin+1], SMeanCH4[zbin+1];
extern double SMeanNH3[zbin+1];
extern double SMeanC2H2[zbin+1], SMeanC2H4[zbin+1], SMeanC2H6[zbin+1], SMeanCH2O2[zbin+1];
extern double SMeanHCN[zbin+1], SMeanN2O[zbin+1], SMeanNO[zbin+1], SMeanNO2[zbin+1], SMeanOCS[zbin+1], SMeanHNO3[zbin+1];
extern int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
extern int    numx, numc, numf, numa, waternum, waterx, numr, numm, numt, nump;
extern double xx[zbin+1][NSP+1];
extern double mkv[], Tnew[], Pnew[];
extern double H2H2CIA[zbin+1][NLAMBDA], H2HeCIA[zbin+1][NLAMBDA], H2HCIA[zbin+1][NLAMBDA], N2H2CIA[zbin+1][NLAMBDA], N2N2CIA[zbin+1][NLAMBDA], CO2CO2CIA[zbin+1][NLAMBDA];
extern double MeanH2H2CIA[], MeanH2HeCIA[], MeanH2HCIA[], MeanN2H2CIA[], MeanN2N2CIA[],MeanCO2CO2CIA[];
extern double SMeanH2H2CIA[], SMeanH2HeCIA[], SMeanH2HCIA[], SMeanN2H2CIA[], SMeanN2N2CIA[], SMeanCO2CO2CIA[];

#endif /* !__GLOBAL_H__ */

/*---- end ------------------------ global.h ---------------------- */
