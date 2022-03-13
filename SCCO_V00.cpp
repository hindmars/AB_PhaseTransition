/*
 * This is the *.cpp file of the Strong Coupling Correction Object (SCCO)
 * 
 * Member functions are declared at SCCO.hpp, and then defined at here.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 


#include <iostream>
#include <cstddef>
#include <cmath>
// using std::istream; using std::ostream;

#include "SCCO_V00.hpp"

//********************************************************************
//********************************************************************
//''' data sheets of  strong coupling corrected material parameters
//'

//constexpr long double SCCO::pi = = 3.14159265358979323846264338328L;
const long double SCCO::c1_arr[18] = {-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413};
const long double SCCO::c2_arr[18] = {-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645};
const long double SCCO::c3_arr[18] = {-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268};
const long double SCCO::c4_arr[18] = {-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518};
const long double SCCO::c5_arr[18] = {-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815};

const long double SCCO::Tc_arr[18] = {0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486}; // mK
const long double SCCO::Ms_arr[18] = {2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82}; // in unit of helium-3 atom
const long double SCCO::VF_arr[18] = {59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23}; // fermi velosity, m.s^-1
const long double SCCO::XI0_arr[18] = {77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76};


//*********************************************************************
//*********************************************************************
//'''     member functions, interfaces of dimensional qualities
//'

long double
SCCO::Tcp(long double p){
  long double Tc = lininterp(Tc_arr, p)*std::pow(10,-3);
  return Tc;
}

long double
SCCO::mEffp(long double p){
  constexpr long double kg = 1.0L;
  const long double u = 1.66053906660L*(std::pow(10.0L,-27))*kg;
  const long double m3 = 3.016293L*u;
  
  long double mEff = lininterp(Ms_arr, p)*m3;;
  return mEff;
}

long double
SCCO::vFp(long double p){
  // unit m.s^-1
  long double vF = lininterp(VF_arr, p);
  return vF;
}

long double
SCCO::xi0p(long double p){
  constexpr long double m = 1.0L;
  const long double nm = (std::pow(10.0L, -9))*m; 
  long double xi0 = lininterp(XI0_arr, p)*nm;
  return xi0;
}  

long double
SCCO::N0p(long double p){
  constexpr long double J = 1.0L, s = 1.0L, pi = 3.14159265358979323846264338328L;;
  const long double hbar = 1.054571817L*(std::pow(10.0L,-34))*J*s;
  long double N0 = (std::pow(mEffp(p),2)*vFp(p))/((2.0L*pi*pi)*std::pow(hbar,3));
  // ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
  return N0;
}


//**********************************************************************
//**********************************************************************
//'''     member functions, interfaces of dimensionless coefficients

long double
SCCO::alpha_bar(long double p, long double T){ return (1.L/3.L)*(T/Tcp(p)-1); }  


long double
SCCO::beta1_bar(long double p, long double T){
  constexpr long double pi = 3.14159265358979323846264338328L;
  const long double zeta3 = std::riemann_zetal(3.0L);
  const long double c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  long double beta1 = c_betai*(-1.0L + (T/Tcp(p))*lininterp(c1_arr, p));

  return beta1;
}  


long double
SCCO::beta2_bar(long double p, long double T){
  constexpr long double pi = 3.14159265358979323846264338328L;
  const long double zeta3 = std::riemann_zetal(3.0L);
  const long double c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  long double beta1 = c_betai*(2.0L + (T/Tcp(p))*lininterp(c2_arr, p));

  return beta1;
}  


long double
SCCO::beta3_bar(long double p, long double T){
  constexpr long double pi = 3.14159265358979323846264338328L;
  const long double zeta3 = std::riemann_zetal(3.0L);
  const long double c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  long double beta1 = c_betai*(2.0L + (T/Tcp(p))*lininterp(c3_arr, p));

  return beta1;
}  


long double
SCCO::beta4_bar(long double p, long double T){
  constexpr long double pi = 3.14159265358979323846264338328L;
  const long double zeta3 = std::riemann_zetal(3.0L);
  const long double c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  long double beta1 = c_betai*(2.0L + (T/Tcp(p))*lininterp(c4_arr, p));

  return beta1;
}


long double
SCCO::beta5_bar(long double p, long double T){
  constexpr long double pi = 3.14159265358979323846264338328L;
  const long double zeta3 = std::riemann_zetal(3.0L);
  const long double c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  long double beta1 = c_betai*(-2.0L + (T/Tcp(p))*lininterp(c5_arr, p));

  return beta1;
}  



//**********************************************************************
//**********************************************************************
//'''               linear intepolation function
//'

long double
SCCO::lininterp(const long double *cX_arr, long double p){
  long double pk, pk1, fp;
  size_t k, k1;

  //  pk = 0.0; pk1 = 2.0; k = 0; k1 = ++k;
  if ((p >= 0.0) && (p < 2.0)) { pk = 0.0L; k = 0; }

  // pk = 2.0; pk1 = 4.0; k = 1; k1 = ++k;
  if ((p >= 2.0) && (p < 4.0)) { pk = 2.0L; k = 1; }
 
  if ((p >= 4.0) && (p < 6.0)) { pk = 4.0L; k = 2; }

  // pk = 6.0; pk1 = 8.0;k = 3; k1 = ++k;
  if ((p >= 6.0) && (p < 8.0)) { pk = 6.0L; k = 3; }

  if ((p >= 8.0) && (p < 10.0)) { pk = 8.0L; k = 4; }

  if ((p >= 10.0) && (p < 12.0)) { pk = 10.0L; k = 5; }

  if ((p >= 12.0) && (p < 14.0)) { pk = 12.0L; k = 6; }

  // pk = 14.0; pk1 = 16.0;k = 7; k1 = ++k;
  if ((p >= 14.0) && (p < 16.0)) { pk = 14.0L; k = 7; }

  if ((p >= 16.0) && (p < 18.0)) { pk = 16.0L; k = 8; }

  if ((p >= 18.0) && (p < 20.0)) { pk = 18.0L; k = 9; }

  if ((p >= 20.0) && (p < 22.0)) { pk = 20.0L; k = 10; }

  if ((p >= 22.0) && (p < 24.0)) { pk = 22.0L; k = 11; }

  if ((p >= 24.0) && (p < 26.0)) { pk = 24.0L; k = 12; }

  if ((p >= 26.0) && (p < 28.0)) { pk = 26.0L; k = 13; }

  if ((p >= 28.0) && (p < 30.0)) { pk = 28.0L; k = 14; }

  if ((p >= 30.0) && (p < 32.0)) { pk = 30.0L; k = 15; }

  if ((p >= 32.0) && (p < 34.0)) { pk = 32.0L; k = 16; }

  fp = ((cX_arr[k+1]-cX_arr[k])/2.0)*(p-pk)+cX_arr[k];
  return fp; 
}









