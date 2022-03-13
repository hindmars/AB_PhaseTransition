/*
 * This is the *.hpp file of the Strong Coupling Correction Object (SCCO)
 * 
 * Member functions are declared at here, and then defined at SCCO.cpp.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 

#ifndef SCCO_H
#define SCCO_H

//#include <string>
#include <iostream>
#include <cstddef>
#include <cmath>


class SCCO {
public:
        SCCO() = default; // constructor

        // interfaces of dimensional qualities
	
        long double Tcp(long double p);
        long double mEffp(long double p);
        long double vFp(long double p);
        long double xi0p(long double p);
        long double N0p(long double p);
  

        // interfaces of dimensionless coeficients; SC-correction parts:
  
        long double alpha_bar(long double p, long double T);
        long double beta1_bar(long double p, long double T);
        long double beta2_bar(long double p, long double T);
        long double beta3_bar(long double p, long double T);
        long double beta4_bar(long double p, long double T);
        long double beta5_bar(long double p, long double T);
        
private:

        // SC-data sheets arries, All associate to SCCO class
        static const long double c1_arr[18]; 
        static const long double c2_arr[18];
        static const long double c3_arr[18];
        static const long double c4_arr[18];
        static const long double c5_arr[18];
        // ***************************************************
        // '''
        static const long double Tc_arr[18];
        static const long double Ms_arr[18];
        static const long double VF_arr[18];
        static const long double XI0_arr[18];

        // linear interpolation member:
        long double lininterp(const long double *cX_arr, long double p);
};

#endif
