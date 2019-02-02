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

#include "params.hpp"
#include "pulsar.hpp"
#include <vector>
#include <random>
#include "create.hpp"
#include "dump.hpp"
#include "evol.hpp"
#include "delete.hpp"
#include "birth.hpp"

int main(int argc, char** argv) {
	//random setup
	std::random_device rd;
	std::mt19937 e2(rd());
	std::uniform_real_distribution<> dist(0, 1);
    
	std::vector<Pulsar> p(Nstart); //main array
    
    dump_init(); //initialize dumping
    
    double I_N = create_all(p, dist, e2); //create initial pulsars
    
    CITable Q_B_citable;
    double I_Q = birth_init(Q_B_citable); //prepare pulsar birth
    
    //calculate timestep from pulsar numbers and distribution function integrals
    double dt = 1e15 * I_N * Nbirth / (M_PI * A * I_Q * Nstart);
    
    for (int i=0; i<=Nsteps; ++i) { //main loop
        if (i % Ndump == 0) dump(p, dt * i, i / Ndump); //do dump every Ndump-th step
        evolve_all(p, dt); //evolve
        delete_all(p); //delete unrelevant pulsars
        birth_all(p, Q_B_citable, dist, e2); //create new pulsar(s)
    }
    
	return 0;
}
