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

//Stuff related to new pulsar generation during the simulation
#include "pulsar.hpp"
#include "params.hpp"
#include <cmath>
#include <vector>
#include <random>
#include "approx.hpp"

//birth functions, not necessarily normed

double Q_P(double P) { return P; }

double Q_chi(double chi) { return M_2_PI; }

double Q_B(double B12) {
	double temp = 1+B12;
	temp *= temp;
	temp *= temp;
	return B12*B12/temp;
}

//generate parameters from random numbers distibuted uniformly in [0,1]

double birth_P(double x) {
	double P2 = Pmin*Pmin*(1-x) + Pmax*Pmax*x;
	return sqrt(P2);
}

double birth_chi(double x) {
	return chimin*(1-x) + chimax*x;
}

double birth_B12(double x, CITable const& Q_B_citable) {
	return invint(x, Q_B_citable);
}

//function that precalculates complicated integrals
double birth_init(CITable& Q_B_citable) {
    Q_B_citable = cumint(Q_B, B12min, B12max, intsteps);
    return Q_B_citable.back().second / 2 * (Pmax*Pmax - Pmin*Pmin);
}

//add new pulsar(s) to the array
void birth_all(std::vector<Pulsar>& p, CITable const& Q_B_citable, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    for (int i=0; i<Nbirth; ++i)
        p.push_back({birth_P(dist(e2)), birth_chi(dist(e2)), birth_B12(dist(e2), Q_B_citable)});
}
