//Approximate operations common for multiple other files
#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <vector>
#include <utility>
#include <functional>

std::vector<std::pair<double, double>> cumint(std::function<double(double)> f, double xmin, double xmax, int nsteps);

double invint(double val, std::vector<std::pair<double, double>> const& citable);

#endif
