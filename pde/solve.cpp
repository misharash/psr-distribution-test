/*
    Copyright 2019 Mykhailo Rashkovetskyi

    This file is part of psr-distribution-test - a program to test
    pulsar distribution solving the kinetic equation as PDE numerically.

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

#include <vector>
#include <cstdio>
#include "params.hpp"
#include "model.hpp"

template<typename T>
int fwrite_vec(std::vector<T> v, FILE* fp) {
    return fwrite(v.data(), sizeof(T), v.size(), fp);
}

void solve_givenB(double B12) {
    //prepare grids
    double dP = (Pmax-Pmin)/NP;
    double dchi = (chimax-chimin)/Nchi;
    std::vector<double> Pgrid(NP+1), chigrid(Nchi+1);
    for (int i=0; i<=NP; ++i)
        Pgrid[i] = Pmin + dP*i;
    for (int i=0; i<=Nchi; ++i)
        chigrid[i] = chimin + dchi*i;
    //prepare file
    system("mkdir dumps");
    char filename[50];
    sprintf(filename, "dumps/dumpB%0.1lf", B12);
    FILE* fp = fopen(filename, "wb");
    //write header
    fprintf(fp, "%d %d\n", NP+1, Nchi+1);
    fwrite(&B12, sizeof(double), 1, fp);
    //write grid coordinates
    fwrite_vec(Pgrid, fp);
    fwrite_vec(chigrid, fp);
    //define two working arrays to do leapfrog
    std::vector<double> warr(Nchi+1, 0), pwarr(Nchi+1, 0);
    //chi-direction flux array
    std::vector<double> fchi(Nchi+1);
    //array of N that will be written after that
    std::vector<double> N(Nchi+1, 0);
    //write two first data rows = zeros
    fwrite_vec(N, fp);
    fwrite_vec(N, fp);
    //main loop, makes a jump from time i-1 to time i+1
    for (int i=1; i<NP; ++i) {
        //fill flux array
        for (int j=1; j<Nchi; ++j)
            fchi[j] = N[j] * chidot(Pgrid[i], chigrid[j], B12);
        //add flux on boundaries
        boundflux(fchi);
        //leapfrog
        for (int j=1; j<Nchi; ++j)
            pwarr[j] += dP * ((fchi[j-1]-fchi[j+1])/dchi + Q(Pgrid[i], chigrid[j], B12));
        //swap arrays without copying
        std::swap(warr, pwarr);
        //calculate N
        for (int j=1; j<Nchi; ++j)
            N[j] = warr[j]/Pdot(Pgrid[i+1], chigrid[j], B12);
        //and write it
        fwrite_vec(N, fp);
    }
    fclose(fp);
}