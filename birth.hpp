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

//Stuff related to new pulsar generation during the simulation
#ifndef BIRTH_HPP
#define BIRTH_HPP

#include "pulsar.hpp"
#include "approx.hpp"
#include <tuple>
#include <random>

std::tuple<double, double, double, double> birth_init(CITable&, CITable&, CITable&, CITable&);

void birth_all(std::vector<Pulsar>&, CITable const&, CITable const&, CITable const&, CITable const&,
                            double, double, double, std::uniform_real_distribution<>&, std::mt19937&);

#endif
