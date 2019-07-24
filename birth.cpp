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
#include "birth.hpp"
#include "params.hpp"
#include <cmath>
#include "approx.hpp"

//birth functions, not necessarily normed

double Q_P(double P) { return 1; }

double Q_chi(double chi) { return sin(chi); }

double Q_B(double B12) { return B12*B12*pow(1+B12, -3.7); }

//`birth` functions on lower P boundary
//cope with lack of pulsar flux through it
//there are two components for usual pulsars

double Q_chibound(double chi) {
    double cchi = cos(chi);
    double schi = sin(chi);
    return (M_PI_2 - chi - schi*cchi)/cchi/cchi/cchi*(1+schi*schi);
}

double Q_Bbound(double B12) { return Q_B(B12); }

//generate parameters from random numbers distibuted uniformly in [0,1]

double birth_P(double x) {
	return Pmin*(1-x) + Pmax*x;
}

double birth_chi(double x) {
    double cchi = cos(chimin)*(1-x)+cos(chimax)*x;
	return acos(cchi);
}

double birth_B12(double x, CITable const& Q_B_citable) {
	return invint(x, Q_B_citable);
}

//the same for boundary `birth`

double birth_chibound(double x, CITable const& Q_chibound1_citable) {
    return invint(x, Q_chibound1_citable);
}

double birth_B12bound(double x, CITable const& Q_B_citable) { //not anything new actually
    return birth_B12(x, Q_B_citable);
}

//function that precalculates complicated integrals
std::tuple<double, double> birth_init(CITable& Q_B_citable, CITable& Q_chibound_citable) {
    Q_B_citable = cumint(Q_B, B12min, B12max, intsteps);
    //Q_Bbound1 is the same as Q_B, skipping
    Q_chibound_citable = cumint(Q_chibound, chimin, chimax, intsteps);
    //Q_Bbound90 is the same as Q_B, skipping
    double I_Q = Q_B_citable.back().second * (Pmax - Pmin);
    double I_Qbound = 0.5 * Pmin * Q_B_citable.back().second * Q_chibound_citable.back().second;
    return std::make_tuple(I_Q, I_Qbound);
}

//add new pulsar(s) to the array
void birth_all(std::vector<Pulsar>& p, CITable const& Q_B_citable, CITable const& Q_chibound_citable,
            double Nbirthb, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    //regular birth
    for (int i=0; i<Nbirth; ++i)
        p.push_back({birth_P(dist(e2)), birth_chi(dist(e2)), birth_B12(dist(e2), Q_B_citable)});
    //boundary `birth`
    if (dist(e2) < Nbirthb)
        p.push_back({Pmin, birth_chibound(dist(e2), Q_chibound_citable), birth_B12bound(dist(e2), Q_B_citable)});
}
