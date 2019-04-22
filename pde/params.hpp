/*
    Copyright 2019 Mykhailo Rashkovetskyi

    This file is part of psr-distribution-test - a program to test
    pulsar distribution solving the kinetic equation as PDE numerically.

    psr-distribution-test is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    psr-distribution-test is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with psr-distribution-test. If not, see
    <https://www.gnu.org/licenses/>.
*/

//Common parameters in one place
#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

//parameters that are assumed to be the same for all pulsars
//const double M = 2e33; //common pulsar mass (gram) - not needed now
const double R = 1e6; //common pulsar radius (cm)

const double c = 3e10; //speed of light

const double A = 1.; //coefficient for Q_BGI (evol.cpp)

const double eps = sqrt(2*M_PI*R/c); //coefficient for C (evol.cpp)

//minimal and maximal periods in seconds
const double Pmin = 0;
const double Pmax = 0.5;
//sections in P direction
const int NP = 100;

//minimal and maximal angles
const double chimin = 0;
const double chimax = M_PI_2;
//sections in chi direction
const int Nchi = 100;

#endif
