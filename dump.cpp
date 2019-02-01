//Data output into file
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "pulsar.hpp"

void dump_init() {
    system("mkdir dumps");
}

void dump(std::vector<Pulsar> const& p, double t, int num) {
    //open file
    char filename[50];
    sprintf(filename, "dumps/dump%05d", num);
    FILE* fp = fopen(filename, "wb");
    
    //write header: time and pulsar count
    fprintf(fp, "%le ", t);
    fprintf(fp, "%d ", p.size());
    fprintf(fp, "\n");
    
    //write all other data in binary format
    for (auto&& psr: p) {
        fwrite(&psr.P, sizeof(double), 1, fp);
        fwrite(&psr.chi, sizeof(double), 1, fp);
        fwrite(&psr.B12, sizeof(double), 1, fp);
    }
}
