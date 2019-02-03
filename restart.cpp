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

//Restart program using data from dump file
#include "pulsar.hpp"
#include <vector>
#include <cstdio>

//function that restores data from dump
void restart(std::vector<Pulsar>& p, int id, double& t, int num, double& dt) {
    //open file
    char filename[50];
    sprintf(filename, "dumps/dump%05d-%02d", num, id);
    FILE* fp = fopen(filename, "rb");
    
    //read header: time, pulsar count and timestep
    char buf[100];
    fgets(buf, 100, fp); //read header string from file to skip to the binary data
    int N; //dump's pulsar count
    sscanf(buf, "%le %d %le", &t, &N, &dt); //parse header string
    
    p.resize(N); //prepare pulsar vector
    //read all other data in binary format
    for (int i=0; i<N; i++) {
        fread(&p[i].P, sizeof(double), 1, fp);
        fread(&p[i].chi, sizeof(double), 1, fp);
        fread(&p[i].B12, sizeof(double), 1, fp);
    }
    
    fclose(fp); //close file
    
    printf("Process %d: restarted from dump %d, time %le s, timestep %le s\n", id, num, t, dt);
}
