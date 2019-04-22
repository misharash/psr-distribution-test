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

//Data output into file
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "pulsar.hpp"

void dump_init() {
    system("mkdir dumps");
}

void dump(std::vector<Pulsar> const& p, int id, double t, int num, double dt) {
    //open file
    char filename[50];
    sprintf(filename, "dumps/dump%05d-%02d", num, id);
    FILE* fp = fopen(filename, "wb");
    
    printf("Process %d: dump %d, time %le s\n", id, num, t);
    
    //write header: time, pulsar count and timestep
    fprintf(fp, "%le ", t);
    fprintf(fp, "%d ", static_cast<int>(p.size()));
    fprintf(fp, "%.16le ", dt); //more precisely since it will be used in calculations after restart
    fprintf(fp, "\n");
    
    //write all other data in binary format
    for (auto&& psr: p) {
        fwrite(&psr.P, sizeof(double), 1, fp);
        fwrite(&psr.chi, sizeof(double), 1, fp);
        fwrite(&psr.B12, sizeof(double), 1, fp);
    }
    
    fclose(fp); //close file
}
