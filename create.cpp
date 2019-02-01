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

//generate parameters from random numbers distributed uniformly in [0,1]

double create_P(double x) {
	double P3 = Pmin*Pmin*Pmin*(1-x) + Pmax*Pmax*Pmax*x;
	return cbrt(P3);
}

double create_chi(double x, std::vector<std::pair<double, double>> const& N_chi_citable) {
    return invint(x, N_chi_citable);
}

double create_B(double x, std::vector<std::pair<double, double>> const& N_B_citable) {
    return invint(x, N_B_citable);
}

//function called from main
void create_all(std::vector<Pulsar>& p, std::uniform_real_distribution<>& dist, std::mt19937& e2) {
    //precalculate complicated integrals
    auto N_chi_citable = cumint(N_chi, chimin, chimax, intsteps);
    auto N_B_citable = cumint(N_B, B12min, B12max, intsteps);
    //generate all parameters for every pulsar
    for (int i=0; i<Nstart; i++) {
        p[i].P = create_P(dist(e2));
        p[i].chi = create_chi(dist(e2), N_chi_citable);
        p[i].B12 = create_B(dist(e2), N_B_citable);
    }
}
