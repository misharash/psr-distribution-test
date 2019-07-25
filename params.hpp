/*
    Copyright 2019 Mykhailo Rashkovetskyi

    This file is part of psr-distribution-test - a program to test
    pulsar distribution using Monte-Carlo method.

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
const double K = 1.; //coefficient for losses (evol.cpp), maybe need to fill it by reasonable number

//simulation settings
const int Nstart = 50000; //how many pulsars to generate initially
const int Nbirth = 1; //how many pulsars to add every birth step
const int Ndump = Nstart/5; //how many steps to do between dumps
const int Nsteps = 50*Nstart; //how many steps to do at all

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
