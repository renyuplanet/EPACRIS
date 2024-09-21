/*----------------------- planet.h --------------------------------

Author: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
Last modified: June 9, 2022
Note: The parameters in this file can be modified to model different planets around different stars
--------------------------------------------------------------------- */

#ifndef _PLANET_H_
#define _PLANET_H_

//--------------------------------------------------------------------- 
#define IN_FILE_NAME    "55cnce-jacob-1.h"
//---------------------------------------------------------------------
#define OUT_DIR         "Results/55cnce/standard_1_test/IW/"
//#define OUT_DIR         "Results/toi270d-valid2-Texp8/"
//#define OUT_DIR         "Results/toi270d-valid3-lowTiso/"
//#define OUT_DIR         "Results/toi270d-valid4-highTiso/"
//#define OUT_DIR         "Results/toi270d-valid5-1DTERRA/"
//#define OUT_DIR         "Results/toi270d-valid6-Tiso2k/"
//---------------------------------------------------------------------
#define COMMENTS        "RC-loop validation"
//#define COMMENTS        "Time Stepping HELIOS, init Guillot; Tint=100K" //to be printed in log file for debugging purposes
//#define COMMENTS        "Dayside: RC test n=10; with RT time stepping N=300; old adiabat and conv scheme" //to be printed in log file for debugging purposes
//#define COMMENTS        "Dayside: RC test n=10; with RT time stepping N=300; New adiabat and New conv scheme" //to be printed in log file for debugging purposes
//--------------------------------------------------------------------- 
//MS, 2022:
#define RadConv_Solver  1   // 0 = Guillot TP profile
                            // 1 = Radiative Convective Climate
//- IF RadConv_Solver = 1 --------------------------------------------- 
#define TIME_STEPPING   0   // 0 = matrix solver for RT fluxes
                            // 1 = time stepping for RT fluxes
#define TS_SCHEME       1   // 0 = testing
                            // 1 = HELIOS scheme (code but not publication)

#define TWO_STR_SOLVER  1   // 0 = 2stream as in Toon+1989 (by Renyu Hu)
                            // 1 = 2 stream as in Heng+2018 (by Markus Scheucher)
#define RT_FLUX_SOLVER  4   // 0 = Thomas Algorithm Bottom-Up on Eq sequence 1 (26,25,26,25,...,26,Bnd.Eq.)
                            // 1 = same as 0 but Top-down
                            // 2 = LU Decomposition (ludcmp&lubksb), Eq sequence 1 as above
                            // 3 = Block-Tridiagonal solver on Eq sequence 2 (TOA,25,26,25,...,25,26,Bnd.Eq.)
                            // 4 = PTRANS-I pentadiagonal solver
                            // 5 = Sogabe 2008 pentadiagonal solver (Algorithms 2&3)

#define CONVEC_ADJUST   1   // 0 = initial dry adiabat (by Renyu Hu)
                            // 1 = 2022 update by Markus Scheucher - see ms_conv_funcs.c
//MS end.
//--------------------------------------------------------------------- 
//--------------------------------------------------------------------- 
/* Initial Concentration Setting */
#define IMODE       0   /*  0: Calculate initial concentrations from chemical equilibrium sub-routines (not rad);
                        1: Import from SPECIES_LIST;
                        2: Import from results of previous calculations
                        3: Calculate initial concentrations from simplied chemical equilibrium formula (not rad);
			4: Import from results of previous calculations in the standard form (TP import only for rad) */

/* Iteration Conditions */
#define NMAX        10       /* Maximum Climate - Chemistry Iterations */
#define NMAX_RC     3       /* Maximum Radiative - Convective Iterations */
#define NMAX_RT     1.0e+3  /* Maximum Radiative Transfer Iterations */
#define NRT_RC      1.0E+3  /*RT steps between Convective adjustments after initial RT equilibrium */

/* Planet Physical Properties */
//M_Jupiter = 1.8982E+27kg (=317.8 M_Earth)
//M_Earth = 5.97237E+24 kg
//------------------------------------
//M_Planet = 8.6 M_E
#define MASS_PLANET 5.134E+25      // kg Planet's mass
//R_Jupiter = 69,911E+3 m (= 10.973*R_Earth)
//R_Earth = 6371E+3 m
//------------------------------------
//R_Planet = 1.95 R_E
#define RADIUS_PLANET 12423450.0     // m Planet's radius
#define ORBIT       0.016       // AU  Planet's semi-major axis, equivalent to around Sun

/* Star spectrum */
#define STAR_SPEC   "Library/Star/55cnc.txt"

/* Atmospheric Properties */
#define RefIdxType  0	/* Type of Refractive Index: 0=Air, 1=CO2, 2=He, 3=N2, 4=NH3, 5=CH4, 6=H2, 7=O2, 8=CO, 9=H2O */
#define AIRM        44.0 /* Initial mean molecular mass of atmosphere, in atomic mass unit*/

#define FADV	    0.6667     /* Advection factor: 0.25=uniformly distributed, 0.6667=no Advection */
//#define THETAREF    0.5000		/* Slant Path Angle in radian */
#define THETAANGLE  30.000		// Slant Path Angle in degrees: 60 for global average, 30 for hemispheric
#define PSURFAB	    0.0000		/* Planet Surface Reflectivity */
#define PSURFEM	    1.0000		/* Planet Surface Emissivity */
/* Initial Heat Flux */
#define TINTSET	    100.0e+0		/* Internal Heat Temperature */


/* Set up calculation grid and temperature-pressure profile intial conditions */
#define TPMODE      0            /* 1: import data from a ZTP list;
                                    0: calculate TP profile from the parametized formula*/
//#define TPLIST      "Library/TPinitial/TP_Jupiter.dat"
//#define TPLIST      "Library/TPinitial/TP_general_exp8.dat"
//#define TPLIST      "Results/toi270d-guillot/NewTemperature.dat"
#define TPLIST      "Results/55cnce/standard_1_test/IW/NewTemperature.dat"
//#define TPLIST      "Results/jupiter-valid1time-Guillot/NewTemperature.dat"

#define zbin        80            /*How many altitude bin?*/

#define PTOP        1.0E-8       /* Pressure at the top of atmosphere in bar */
#define TTOP        2100.0       /* Temperature at PTOP in K */
#define TINV        0            /* set to 1 if there is a temperature inversion */
#define PMIDDLE     1.0E-3       /* Pressure at the top of troposphere or stratosphere (if inversion) */
#define TMIDDLE     2100.0       /* Tempeature at PMIDDLE */
#define PSTR        1.0E-1       /* Pressure at the top of troposphere (only used if inversion) */
#define TSTR        2100.0       /* Temperature at PSTR */
#define PTROP       2.0E-1       /* Pressure at the bottom of troposphere */
#define TTROP       2100.0       /* Temperature at PTROP */
#define PBOTTOM     1.0E+0       /* Pressure at the bottom of atmosphere */

/* Molecular Species */
#define SPECIES_LIST "Condition/SpeciesList/species_general_CHONS_exp.dat"
//- CONVECTION -----------------------
#define NCONDENSIBLES 2  //how many potentially condensing species for moist adiabat, cloud formation, and rain-out
#define CONDENSIBLES (int[]){7,9} //H2O=7; NH3=9; CO=20; CH4=21; CO2=52; H2=53; O2=54; N2=55

/* Reaction List */
#define REACTION_LIST "Condition/ReactionList/zone_general_CHO.dat"

/* Radiative Convective Calculation*/
#define IFRC        1        /* do we update radiative-convective boundary in each radiative balance iteration? */

#define Tol_RC_R    1.0E-2  /* convergence tolerance in unbalanced radiative flux, per internal heat flux (net outgoing flux) */
#define Tol_RC      1.0E+0 /* convergence tolerance in unbalanced radiative flux, in absolute quantity in W/m2 (satisfying any of them is ok) */
#define Tol_FRATIO  1.0E-7  //convergence tolerance in unbalanced radiative flux per layer against layer radiance SIGMA*T^4

#define R_RELAX     1.0e-1     /* relaxation factor in implicit Euler stepping of radiative balance */
#define DT_MAX      1.0e-1     /* maximum dT change in implicit Euler stepping of radiative balance */
#define RJACOB_TEMPVAR 1.0e+0  /* temperature variation in calculating the Jacobian */

/* Input choices for the infrared opacities */
/* Must be set to the same as the opacity code */

#define CROSSHEADING		"Library/Opacity/H2_FullT_LowRes/"

#define NTEMP       20             /* Number of temperature points in grid   */
#define TLOW        100.0           /* Temperature range in K                 */
#define THIGH       2000.0

#define NPRESSURE   10         /* Number of pressure points in grid      */
#define PLOW        1.0e-01         /* Pressure range in Pa                   */
#define PHIGH       1.0e+08

#define NLAMBDA     16000         /* Number of wavelength points in grid    */
#define LAMBDALOW   1.0e-07    /* Wavelength range in m                  */
#define LAMBDAHIGH  2.0e-04
#define LAMBDATYPE  1         /* LAMBDATYPE=1 -> constant resolution    */
                             /* LAMBDATYPE=2 -> constant wave step     */


/* Parameters to be adapted in overal model upgrade */
#define NSP         111                 /*Number of species in the standard list*/
#define NKin        645   /*Number of Regular Chemical Reaction in the standard list*/
#define NKinM       90  /*Number of Thermolecular Reaction in the standard list*/
#define NKinT       93  /*Number of Thermal Dissociation Reaction in the standard list*/
#define NPho        71   /*Number of Photochemical Reaction in the standard list*/
#define	THREEBODY   1.0	/* Enhancement of THREEBODY Reaction when CO2 dominant */
#define NATOMS      23            /* Number of atoms for chemical equil     */
#define NMOLECULES  172       /* Number of molecules for chemical equil */
#define MOL_DATA_FILE "Library/ChemEqu/molecules_all.dat" /* Data file for chemical equilibrium calculation */


/* All following parameters are rarely used, set to the default values */
#define STAR_TEMP   157.18//GREY ONLY     /* Irradiance Temperature at 1 AU */
#define FaintSun    1.0				/* Faint early Sun factor */
#define OUT_FILE1   "AuxillaryOut/Conx.dat"
#define OUT_FILE2   "AuxillaryOut/Conf.dat"
#define IFIMPORTH2O 0		/* When H2O is set to constant, 1=import mixing ratios */
#define IFIMPORTCO2 0		/* When H2O is set to constant, 1=import mixing ratios */
#define AERRADFILE1 "Library/Aerosol/H2SO4AER.dat"	/* radiative properties of H2SO4 */
#define AERRADFILE2 "Library/Aerosol/S8AER.dat"	/* radiative properties of S8 */

/* Planet Orbital Properties */
#define TIDELOCK    0 //GREY ONLY   	/* If the planet is tidelly locked */
#define PAB         0.343 //GREY ONLY	/* Planet Bond Albedo */
#define DELADJUST   1			/* Whether use the delta adjustment in the 2-stream diffuse radiation */
#define TAUTHRESHOLD 0.1			/* Optical Depth Threshold for multi-layer diffuse radiation */
#define TAUMAX	    1000.0		/* Maximum optical Depth in the diffuse radiation */
#define TAUMAX1	    1000.0		/* Maximum optical Depth in the diffuse radiation */
#define TAUMAX2	    1000000.0
#define IFDIFFUSE   1			/* Set to 1 if want to include diffuse solar radiation into the photolysis rate */

#define IFUVMULT    0			/* Whether do the UV Multiplying */
#define FUVMULT	    1.0E+3		/* Multiplying factor for FUV radiation <200 nm */
#define MUVMULT	    1.0E+2		/* Multiplying factor for MUV radiation 200 - 300 nm */
#define NUVMULT	    1.0E+1		/* Multiplying factor for NUV radiation 300 - 400 nm */


#define WaveMax     1000.0 /*Maximum Wavelength in nm for the Calculation of UV-visible radiation and photolysis rates*/
#define TDEPMAX	    300.0 /* Maximum Temperature-dependence Validity for UV Cross sections */
#define TDEPMIN     200.0 /* Minimum Temperature-dependence Validity for UV Cross sections */

/* The criteria of convergence */
#define Tol1        1.0E+10
#define Tol2        1.0E-16

/* Mode of iteration */
#define	TSINI	    1.0E-18	/* Initial Trial Timestep, generally 1.0E-8 */
#define FINE1       1	/* Set to one for fine iteration: Set to 2 to disregard the bottom boundary layers */
#define FINE2       1 /* Set to one for fine iteration: Set to 2 to disregard the fastest varying point */
#define TMAX        1.0E+12 /* Maximum of time step */
#define TMIN        1.0E-5 /* Minimum of time step */
#define TSPEED	    1.0E+12 /* Speed up factor */

#define NMAXT	    1.0E+18 /* Maximum iteration cumulative time in seconds */
#define MINNUM      1.0E-0 /* Minimum number density in denominator */

/* Aerosol Species */
#define AIRVIS	    1.0E-5	/* Dynamic viscosity in SI*/
#define	AERSIZE	    1.0E-7	/* diameter in m */
#define AERDEN	    1.84E+3	/* density in SI */
#define	NCONDEN	    1	/* Calculate the condensation every NCONDEN iterations */
#define IFGREYAER   0	/* Contribute to the grey atmosphere Temperature? 0=no, 1=yes */
#define SATURATIONREDUCTION 0.2 /* Ad hoc reduction factor for saturation pressure of water */



/* Parametization of Eddy Diffusion Coefficient */
#define EDDYPARA    1	/* =1 from Parametization, =2 from imported list */
#define KET         1.0E+6 /*unit cm2 s-1*/
#define KEH         1.0E+6
#define ZT          200.0  /*unit km*/
#define Tback       1E+4
#define KET1        1.0E+6
#define KEH1        1.0E+8
#define EDDYIMPORT  "Data/EddyH2.dat"
#define MDIFF_H_1   4.87
#define MDIFF_H_2   0.698
#define MDIFF_H2_1  2.80
#define MDIFF_H2_2  0.740
#define MDIFF_H2_F  1.0

/* Parameters of rainout rates */
#define	RainF	    0.0	/* Rainout factor, 0 for no rainout, 1 for earthlike normal rainout, <1 for reduced rainout */
#define	CloudDen    1.0	/* Cloud density in the unit of g m-3 */

/* Output Options */
#define NPRINT      1E+2               /* Printout results and histories every NPRINT iterations */
#define HISTORYPRINT 0			/* print out time series of chemical composition if set to 1 */


/* IR emission spectra output options */

#define IRLamMin    1.0		/* Minimum wavelength in the IR emission output, in microns */
#define IRLamMax    100.0	/* Maximum wavelength in the IR emission output, in microns */
#define IRLamBin    9999		/* Number of wavelength bin in the IR emission spectra */
#define Var1STD	    43
#define	Var2STD	    45
#define Var3STD	    78
#define	Var4STD	    111
#define Var1RATIO   0.0
#define	Var2RATIO   0.0
#define Var3RATIO   0.0
#define	Var4RATIO   0.0

/*  Stellar Light Reflection output options */
#define UVRFILE	    "Result/GJ1214_Figure/Reflection.dat"  /* Output spectrum file name */
#define UVRFILEVar1 "Result/GJ1214_Figure/ReflectionVar1.dat"  /* Output spectrum file name */
#define UVRFILEVar2 "Result/GJ1214_Figure/ReflectionVar2.dat"  /* Output spectrum file name */
#define UVRFILEVar3 "Result/GJ1214_Figure/ReflectionVar3.dat"  /* Output spectrum file name */
#define UVRFILEVar4 "Result/GJ1214_Figure/ReflectionVar4.dat"  /* Output spectrum file name */
#define UVROPTFILE  "Result/GJ1214_Figure/UVROpt.dat" /* Output spectrum file name*/

/* Stellar Light Transmission output options */
#define UVTFILE	    "Result/GJ1214_Figure/Transmission.dat" /* Output spectrum file name */
#define UVTFILEVar1 "Result/GJ1214_Figure/TransmissionVar1.dat" /* Output spectrum file name */
#define UVTFILEVar2 "Result/GJ1214_Figure/TransmissionVar2.dat" /* Output spectrum file name */
#define UVTFILEVar3 "Result/GJ1214_Figure/TransmissionVar3.dat" /* Output spectrum file name */
#define UVTFILEVar4 "Result/GJ1214_Figure/TransmissionVar4.dat" /* Output spectrum file name */
#define UVTOPTFILE  "Result/GJ1214_Figure/UVTOpt.dat" /* Output spectrum file name*/

/* Thermal Emission output options */
#define IRFILE	    "Result/GJ1214_Figure/Emission.dat"	     /* Output spectrum file name */
#define IRFILEVar1  "Result/GJ1214_Figure/EmissionVar1.dat"	 /* Output spectrum file name */
#define IRFILEVar2  "Result/GJ1214_Figure/EmissionVar2.dat"	 /* Output spectrum file name */
#define IRFILEVar3  "Result/GJ1214_Figure/EmissionVar3.dat"	 /* Output spectrum file name */
#define IRFILEVar4  "Result/GJ1214_Figure/EmissionVar4.dat"	 /* Output spectrum file name */
#define IRCLOUDFILE "Result/GJ1214_Figure/CloudTopE.dat"      /* Output emission cloud top file name */

/* Cloud Top Determination */
#define OptCloudTop 1.0	/* Optical Depth of the Cloud Top */

#endif

/* 1 Tg yr-1 = 3.7257E+9 H /cm2/s for earth */
