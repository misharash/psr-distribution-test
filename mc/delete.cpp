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

//Stuff related to removing pulsars
#include "pulsar.hpp"
#include "params.hpp"
#include <vector>
#include <algorithm>

//checks if the pulsar crossed parameter boundaries
bool bound_check(Pulsar& psr) {
    //chi boundary check - their crossing is unphysical so just fix them
    if (psr.chi > chimax) psr.chi = chimax;
    if (psr.chi < chimin) psr.chi = chimin;
    //P boundary check
    return (psr.P > Pmax)||(psr.P < Pmin); //second part shouldn't apply in any braking model
}

//function that deletes pulsars from array
void delete_all(std::vector<Pulsar>& p) {
    p.erase(std::remove_if(p.begin(), p.end(), bound_check), p.end()); //deletes pulsar that crossed boundaries
    //maybe later should apply here death function
}
