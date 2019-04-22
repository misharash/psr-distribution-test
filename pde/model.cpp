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

//Functions that depend on chosen model
#include "params.hpp"
#include <cmath>
#include <vector>

//Period time-derivative, times 1e-15 omitted
double Pdot(double P, double chi, double B12) {
    return A*pow(B12, 10./7)*pow(P, 1./14)*pow(cos(chi), 1.5)
    + eps*B12*B12/pow(P, 1.5);
}

//Inclination angle time-derivative, times 1e-15/s omitted
double chidot(double P, double chi, double B12) {
    return A*pow(B12, 10./7)*sqrt(cos(chi))*sin(chi)/pow(P, 13./14);
}

//set flux on chi boundaries
void boundflux(std::vector<double> fchi) {
    fchi[0] = 0;
    fchi[fchi.size()-1] = fchi[fchi.size()-2];
}

//birth functions, not normed since it can be done later
double Q_P(double P) { return P; }

double Q_chi(double chi) { return 2; }

double Q_B(double B12) {
	double temp = 1+B12;
	temp *= temp;
	temp *= temp;
	return B12*B12/temp;
}

double Q(double P, double chi, double B12) {
    return Q_P(P)*Q_chi(chi)*Q_B(B12);
}
