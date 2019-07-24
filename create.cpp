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

//Stuff related to pulsar creation in the beginning of the simulation
#include "pulsar.hpp"
#include "params.hpp"
#include "approx.hpp"
#include <random>

//initial distribution functions, not neccesarily normed

double N_P(double P) { return P*P; }

double N_chi(double chi) {
	double schi = sin(chi);
	double cchi = cos(chi);
	return (M_PI_2-chi-schi*cchi)/cchi/cchi/cchi;
}

double N_B(double B12) { return pow(1+B12, -3.7); }

//generate parameters from random numbers distributed uniformly in [0,1]

double create_P(double x) {
	double P3 = Pmin*Pmin*Pmin*(1-x) + Pmax*Pmax*Pmax*x;
	return cbrt(P3);
}

double create_chi(double x, CITable& N_chi_citable) {
    return invint(x, N_chi_citable);
}

double create_B(double x) {
    double Bp = (1-x)*pow(1+B12min, -2.7) + x*pow(1+B12max, -2.7);
    return pow(Bp, -1/2.7)-1;
}

//function that creates initial pulsars and returns distribution function integral over parameter area
double create_all(std::vector<Pulsar>& p, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    //precalculate complicated integrals
    auto N_chi_citable = cumint(N_chi, chimin, chimax, intsteps);
    //regular pulsar count
    int Nreg = Nstart;
    //generate all parameters for every regular pulsar
    for (int i=0; i<Nreg; ++i) {
        p[i].P = create_P(dist(e2));
        p[i].chi = create_chi(dist(e2), N_chi_citable);
        p[i].B12 = create_B(dist(e2));
    }
    return N_chi_citable.back().second * (Pmax*Pmax*Pmax - Pmin*Pmin*Pmin)/3 * (pow(1+B12min, -2.7) - pow(1+B12max, -2.7))/2.7;
}
