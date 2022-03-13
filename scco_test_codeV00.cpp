
#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "SCCO_V00.hpp"


int main(){

  long double T = 1.5*std::pow(10,-3);
  long double p = 32.0;
  long double something = 3.14159265358979323846264338328L;

  SCCO SC;

  std::cout << SC.beta1_bar(p, T) << " and Tc is " << SC.Tcp(p) << std::endl;
  std::cout << something << SC.alpha_bar(p, T) << std::endl;
  std::cout << SC.N0p(p) << std::endl;

  return 0;


}  
