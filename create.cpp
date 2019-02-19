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
	return (chi - schi*cchi)/(schi*schi*schi*sqrt(cchi));
}

double N_B(double B12) {
	double temp = 1+B12;
	temp *= temp;
	temp *= temp;
	return pow(B12, 4./7)/temp;
}

//initial distribution functions for orthogonal, not neccesarily normed

double N_P90(double P) { return pow(P, 25./7); }

double N_B90(double B12) {
	double temp = 1+B12;
	temp *= temp;
	temp *= temp;
	return 1/temp;
}

//generate parameters from random numbers distributed uniformly in [0,1]

double create_P(double x) {
	double P3 = Pmin*Pmin*Pmin*(1-x) + Pmax*Pmax*Pmax*x;
	return cbrt(P3);
}

double create_chi(double x, CITable& N_chi_citable) {
    return invint(x, N_chi_citable);
}

double create_B(double x, CITable& N_B_citable) {
    return invint(x, N_B_citable);
}

double create_P90(double x) {
    double Pp = pow(Pmin, 32./7)*(1-x) + pow(Pmax, 32./7)*x;
    return pow(Pp, 7./32);
}

double create_B90(double x) {
    double temp1 = 1+B12min;
    double temp2 = 1+B12max;
    double Bp = (1-x)/temp1/temp1/temp1 + x/temp2/temp2/temp2;
    return 1/cbrt(Bp)-1;
}

//function that creates initial pulsars and returns distribution function integral over parameter area
double create_all(std::vector<Pulsar>& p, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    //precalculate complicated integrals
    auto N_chi_citable = cumint(N_chi, chimin, chimax, intsteps);
    auto N_B_citable = cumint(N_B, B12min, B12max, intsteps);
    double norm = N_chi_citable.back().second * N_B_citable.back().second / 3 * (Pmax*Pmax*Pmax - Pmin*Pmin*Pmin);
    //similar integrals for orthogonal pulsars
    double temp1 = 1+B12min;
    double temp2 = 1+B12max;
    double norm90 = 7.*M_PI/29/eps * A * 7./32 * (pow(Pmax, 32./7) - pow(Pmin, 32./7)) / 3 * (1/temp1/temp1/temp1 - 1/temp2/temp2/temp2);
    //regular pulsar count
    int Nreg = static_cast<int>(norm / (norm + norm90) * Nstart);
    //generate all parameters for every regular pulsar
    for (int i=0; i<Nreg; ++i) {
        p[i].P = create_P(dist(e2));
        p[i].chi = create_chi(dist(e2), N_chi_citable);
        p[i].B12 = create_B(dist(e2), N_B_citable);
    }
    //add orthogonal pulsars
    for (int i=Nreg; i<Nstart; ++i) {
        p[i].P = create_P90(dist(e2));
        p[i].chi = chimax;
        p[i].B12 = create_B90(dist(e2));
    }
    return norm + norm90;
}
