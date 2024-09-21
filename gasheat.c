/* Function to calculate the heat capacity of gas as a function of temperature */
/* output unit J MOLE-1 K-1 */
/* input unit: K */

#include <math.h>

double AirHeat(double ); /* Type 0 */
double CO2Heat(double ); /* Type 1 */
double HeHeat(double );  /* Type 2 */
double N2Heat(double );  /* Type 3 */
double NH3Heat(double ); /* Type 4 */
double CH4Heat(double ); /* Type 5 */
double H2Heat(double );  /* Type 6 */
double O2Heat(double );  /* Type 7 */
double COHeat(double );  /* Type 8 */
double H2OHeat(double );  /* Type 9 */

double AirHeat(double T)
{
	double cp1, cp2, cp;
	cp1 = N2Heat(T);
	cp2 = O2Heat(T);
	cp  = cp1*0.8 + cp2*0.2;
	return cp;
}

double CO2Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<298.0) {
		T1=298.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	
	t=T1/1000.0;
	if (T1<1200.0) {
		A=24.99735;
		B=55.18696;
		C=-33.69137;
		D=7.948387;
		E=-0.136638;
	} else {
		A=58.16639;
		B=2.720074;
		C=-0.492289;
		D=0.038844;
		E=-6.447293;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double HeHeat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<298.0) {
		T1=298.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	t=T1/1000.0;
	A=20.78603;
	B=4.850638E-10;
	C=-1.582916E-10;
	D=1.525102E-11;
	E=3.196347E-11;
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double N2Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<100.0) {
		T1=100.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	
	t=T1/1000.0;
	if (T1<500.0) {
		A=28.98641;
		B=1.853978;
		C=-9.647459;
		D=16.63537;
		E=0.000117;
	} else if (T1<2000.0) {
		A=19.50583;
		B=19.88705;
		C=-8.598535;
		D=1.369784;
		E=0.527601;
	} else {
		A=35.51872;
		B=1.128728;
		C=-0.196103;
		D=0.014662;
		E=-4.553760;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double NH3Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<298.0) {
		T1=298.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	
	t=T1/1000.0;
	if (T1<1400.0) {
		A=19.99563;
		B=49.77119;
		C=-15.37599;
		D=1.921168;
		E=0.189174;
	} else {
		A=52.02427;
		B=18.48801;
		C=-3.765128;
		D=0.248541;
		E=-12.45799;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double CH4Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<298.0) {
		T1=298.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	
	t=T1/1000.0;
	if (T1<1300.0) {
		A=-0.703029;
		B=108.4773;
		C=-42.52157;
		D=5.862788;
		E=0.678565;
	} else {
		A=85.81217;
		B=11.26467;
		C=-2.114146;
		D=0.138190;
		E=-26.42221;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double H2Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<298.0) {
		T1=298.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}

	t=T1/1000.0;
	if (T1<1000.0) {
		A=33.066178;
		B=-11.363417;
		C=11.432816;
		D=-2.772874;
		E=-0.158558;
	} else if (T1<2500.0) {
		A=18.563083;
		B=12.257357;
		C=-2.859786;
		D=0.268238;
		E=1.97799;
	} else {
		A=43.41356;
		B=-4.293079;
		C=1.272428;
		D=-0.096876;
		E=-20.533862;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}

double O2Heat(double T)
{
	double t, A, B, C, D, E, cp;
	double T1;
	if (T<100.0) {
		T1=100.0;
	} else if (T>6000.0) {
		T1=6000.0;
	} else {
		T1=T;
	}
	
	t=T1/1000.0;
	if (T1<700.0) {
		A=31.32234;
		B=-20.23531;
		C=57.86644;
		D=-36.50624;
		E=-0.007374;
	} else if (T1<2000.0) {
		A=30.03235;
		B=8.772972;
		C=-3.988133;
		D=0.788313;
		E=-0.741599;
	} else {
		A=20.91111;
		B=10.72071;
		C=-2.020498;
		D=0.146449;
		E=9.245722;
	}
	cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
	return cp;
}


double COHeat(double T)
{
    double t, A, B, C, D, E, cp;
    double T1;
    if (T<298.0) {
        T1=298.0;
    } else if (T>6000.0) {
        T1=6000.0;
    } else {
        T1=T;
    }
    
    t=T1/1000.0;
    if (T1<1300.0) {
        A=25.56759;
        B=6.096130;
        C=4.054656;
        D=-2.671301;
        E=0.131021;
    } else {
        A=35.15070;
        B=1.300095;
        C=-0.205921;
        D=0.013550;
        E=-3.282780;
    }
    cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
    return cp;
}


double H2OHeat(double T)
{
    double t, A, B, C, D, E, cp;
    double T1;
    if (T<500.0) {
        T1=500.0;
    } else if (T>6000.0) {
        T1=6000.0;
    } else {
        T1=T;
    }
    
    t=T1/1000.0;
    if (T1<1700.0) {
        A=30.09200;
        B=6.832514;
        C=6.793435;
        D=-2.534480;
        E=0.082139;
    } else {
        A=41.96426;
        B=8.622053;
        C=-1.499780;
        D=0.098119;
        E=-11.15764;
    }
    cp = A + B*t +  C*t*t + D*t*t*t + E/t/t;
    return cp;
}
