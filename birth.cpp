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

double Q_P(double P) { return P; }

double Q_chi(double chi) { return 2; }

double Q_B(double B12) {
	double temp = 1+B12;
	temp *= temp;
	temp *= temp;
	return B12*B12/temp;
}

//`birth` functions on lower P boundary
//cope with lack of pulsar flux through it
//there are two components for usual pulsars

double Q_chibound1(double chi) {
    double cchi = cos(chi);
    double schi = sin(chi);
    return (chi - schi * cchi)*cchi/schi/schi/schi;
}

double Q_Bbound1(double B12) { return Q_B(B12); }

double Q_chibound2(double chi) {
    double cchi = cos(chi);
    double schi = sin(chi);
    return (chi - schi * cchi)/sqrt(cchi)/schi/schi/schi;
}

double Q_Bbound2(double B12) { return pow(B12, 4./7) * Q_B(B12); }

//separate function for orthogonal, though not actually new

double Q_Bbound90(double B12) { return Q_B(B12); }

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

//the same for boundary `birth`

double birth_chibound1(double x, CITable const& Q_chibound1_citable) {
    return invint(x, Q_chibound1_citable);
}

double birth_B12bound1(double x, CITable const& Q_B_citable) { //not anything new actually
    return birth_B12(x, Q_B_citable);
}

double birth_chibound2(double x, CITable const& Q_chibound2_citable) {
    return invint(x, Q_chibound2_citable);
}

double birth_B12bound2(double x, CITable const& Q_Bbound2_citable) {
    return invint(x, Q_Bbound2_citable);
}

double birth_B12bound90(double x, CITable const& Q_B_citable) { //not anything new actually
    return birth_B12(x, Q_B_citable);
}

//function that precalculates complicated integrals
std::tuple<double, double, double, double> birth_init(CITable& Q_B_citable, CITable& Q_chibound1_citable,
                    CITable& Q_chibound2_citable, CITable& Q_Bbound2_citable) {
    Q_B_citable = cumint(Q_B, B12min, B12max, intsteps);
    //Q_Bbound1 is the same as Q_B, skipping
    Q_Bbound2_citable = cumint(Q_Bbound2, B12min, B12max, intsteps);
    Q_chibound1_citable = cumint(Q_chibound1, chimin, chimax, intsteps);
    Q_chibound2_citable = cumint(Q_chibound2, chimin, chimax, intsteps);
    //Q_Bbound90 is the same as Q_B, skipping
    double I_Q = 1e-15 * M_PI * A * Q_B_citable.back().second / 2 * (Pmax*Pmax - Pmin*Pmin);
    double I_Qbound1 = 1e-15 * A * pow(Pmin, 29./14) * Q_B_citable.back().second * Q_chibound1_citable.back().second;
    double I_Qbound2 = 1e-15 * eps * sqrt(Pmin) * Q_Bbound2_citable.back().second * Q_chibound2_citable.back().second;
    double I_Qbound90 = 1e-15 * 7*M_PI/29 * pow(Pmin, 37./14) * Q_B_citable.back().second;
    return std::make_tuple(I_Q, I_Qbound1, I_Qbound2, I_Qbound90);
}

//add new pulsar(s) to the array
void birth_all(std::vector<Pulsar>& p, CITable const& Q_B_citable, CITable const& Q_chibound1_citable,
            CITable const& Q_chibound2_citable, CITable const& Q_Bbound2_citable, double Nbirthb1,
            double Nbirthb2, double Nbirthb90, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    //regular birth
    for (int i=0; i<Nbirth; ++i)
        p.push_back({birth_P(dist(e2)), birth_chi(dist(e2)), birth_B12(dist(e2), Q_B_citable)});
    //boundary `birth`
    if (dist(e2) < Nbirthb1)
        p.push_back({Pmin, birth_chibound1(dist(e2), Q_chibound1_citable), birth_B12bound1(dist(e2), Q_B_citable)});
    if (dist(e2) < Nbirthb2)
        p.push_back({Pmin, birth_chibound2(dist(e2), Q_chibound2_citable), birth_B12bound1(dist(e2), Q_Bbound2_citable)});
    if (dist(e2) < Nbirthb90)
        p.push_back({Pmin, chimax, birth_B12bound90(dist(e2), Q_B_citable)});
}
