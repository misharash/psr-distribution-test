# Copyright 2019 Mykhailo Rashkovetskyi
#
# This file is part of psr-distribution-test - a program to test
# pulsar distribution using Monte-Carlo method.
#
# psr-distribution-test is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# psr-distribution-test is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with psr-distribution-test. If not, see
# <https://www.gnu.org/licenses/>.

from numpy import fromfile, float64


# function that reads data from single dump
def rds(filename):
    f = open(filename, "rb")
    # read header
    header = f.readline().decode()
    L = header.split()
    NP = int(L[0])-1
    Nchi = int(L[1])-1
    # read body
    body = fromfile(f, dtype=float64, count=(NP+2)*(Nchi+2))
    B12, body = body[0], body[1:]
    P, body = body[:NP+1], body[NP+1:]
    chi, body = body[:Nchi+1], body[Nchi+1:]
    N = body.reshape(NP+1, Nchi+1)
    return P, chi, B12, N
