//Common parameters in one place
#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

//parameters that are assumed to be the same for all pulsars
const double M = 2e33; //common pulsar mass (gram)
const double R = 1e6; //common pulsar radius (cm)

const double c = 3e10; //speed of light

//simulation settings
const int Nstart = 1000000; //how many pulsars to generate initially
const int Nbirth = 1; //how many pulsars to add every birth step

//numerical settings
const int intsteps = 100000; //how many steps to use in numerical integration

//minimal and maximal periods in seconds
const double Pmin = 0.03;
const double Pmax = 0.5;

//minimal and maximal angles
const double chimin = 0;
const double chimax = M_PI_2;

//minimal and maximal magnetic fields in 10^12 G
const double B12min = 1e-2;
const double B12max = 10;

#endif
