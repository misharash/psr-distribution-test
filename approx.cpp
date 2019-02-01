//Approximate operations common for multiple other files
#include <vector>
#include <utility>
#include <functional>
#include <algorithm>

//dumb middle-rectangles cumulative integrator
std::vector<std::pair<double, double>> cumint(std::function<double(double)> f, double xmin, double xmax, int nsteps) {
	std::vector<std::pair<double, double>> r(nsteps+1);
	r[0].first = xmin;
	r[0].second = 0;
	double dx = (xmax - xmin) / nsteps;
	for (int i=1; i<=nsteps; i++) {
		r[i].first = xmin + dx * i;
		r[i].second = r[i-1].second + dx * f(r[i].first - dx/2);
	}
	r[nsteps].first = xmax;
	return r;
}

//dumb function that inverses normed integral using cumulative table from above function
double invint(double val, std::vector<std::pair<double, double>> const& citable) {
    val *= citable.back().second; //applying norm here
    //binary search to find first element with integral more or equal than value
    auto nextIt = std::lower_bound(citable.begin(), citable.end(), val, [](std::pair<double, double> const& x, double y) { return x.second < y; });
    //boundary checks that I hope won't apply
    if (nextIt == citable.begin()) ++nextIt;
    if (nextIt == citable.end()) --nextIt;
    auto prevIt = nextIt - 1;
    double frac = (val - prevIt->second) / (nextIt->second - prevIt->second);
    double inv = prevIt->first + frac * (nextIt->first - prevIt->first);
    //boundary checks that I hope won't apply
    if (inv < citable.front().first) inv = citable.front().first;
    if (inv > citable.back().first) inv = citable.back().first;
    return inv;
}
