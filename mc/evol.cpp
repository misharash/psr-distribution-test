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

//Pulsar evolution with time
#include "pulsar.hpp"
#include "params.hpp"
#include <cmath>
#include <vector>

void evolve(Pulsar& p, double dt) {
	double schi = sin(p.chi);
	double cchi = cos(p.chi);
	double Q_BGI = 0;
	if (cchi > 0) //prevent division by zero or square root of negative
		Q_BGI = A * pow(p.P, 15./14) * pow(p.B12, -4./7) / sqrt(cchi);
	double C = eps/sqrt(p.P);
	p.P += dt*1e-15 * p.B12*p.B12/p.P * (Q_BGI * cchi*cchi + C);
	p.chi += dt*1e-15 * Q_BGI * p.B12*p.B12/p.P/p.P * schi * cchi;
}

void evolve_all(std::vector<Pulsar>& p, double dt) {
    for (auto&& psr: p)
        evolve(psr, dt);
}
