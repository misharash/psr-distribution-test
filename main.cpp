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
#include "restart.hpp"
#include "dump.hpp"
#include "evol.hpp"
#include "delete.hpp"
#include "birth.hpp"
#include <cfenv> //nan handling
#include <unistd.h> //forking
#include <sys/wait.h> //waiting for children
#include <cstdio>

//first argument is number of processes to create
//second argument if present is number of dump to try to restart from
int main(int argc, char** argv) {
    //stuff common for all processes
    feenableexcept(FE_INVALID | FE_OVERFLOW); //raise exception when nan creates
    
    CITable Q_B_citable;
    double I_Q = birth_init(Q_B_citable); //prepare pulsar birth
    
    dump_init(); //initialize dumping
    
    //get number of processes to create, default is 1
    int nproc = 1;
    if (argc > 1) sscanf(argv[1], "%d", &nproc);
    
    int myid; //process individual id
    for (myid = nproc-1; myid > 0; --myid) //main process myid will be 0
        if (!fork()) break; //only main process creates children
    printf("Process with id %d started\n", myid);
    
    //stuff individual for every process
	//random setup
	std::random_device rd;
	std::mt19937 e2(rd());
	std::uniform_real_distribution<> dist(0, 1);
    
	std::vector<Pulsar> p(Nstart); //main array
    double dt; //timestep
    
    //restart-related stuff
    int nrdump = 0; //number of first dump/dump to restart from
    double t0 = 0; //time of start/restart
    if (argc > 2) { //restart only if second argument exists
        sscanf(argv[2], "%d", &nrdump);
        restart(p, myid, t0, nrdump, dt); //attempt to read needed data from dump
    }
    else {
        double I_N = create_all(p, dist, e2); //create initial pulsars
        //calculate timestep from pulsar numbers and distribution function integrals
        dt = 1e15 * I_N * Nbirth / (M_PI * A * I_Q * Nstart);
        printf("Process %d: timestep %le s\n", myid, dt);
    }
    
    for (int i=0; i<=Nsteps; ++i) { //main loop
        if (i % Ndump == 0) dump(p, myid, t0 + dt * i, nrdump + i / Ndump, dt); //do dump every Ndump-th step
        evolve_all(p, dt); //evolve
        delete_all(p); //delete unrelevant pulsars
        birth_all(p, Q_B_citable, dist, e2); //create new pulsar(s)
    }
    
    int wstatus; //needed to call wait
    if (!myid) //main process waits for all children
        for (int i=1; i<nproc; ++i) wait(&wstatus);
    
    printf("Process with id %d finished\n", myid);
	return 0;
}
