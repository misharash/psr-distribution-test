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

from numpy import fromfile, float64, zeros, append
from glob import glob


# function that reads data from single dump
def rds(filename):
    f = open(filename, "rb")
    # read header
    header = f.readline().decode()
    L = header.split()
    t = float(L[0])
    N = int(L[1])
    # read body
    body = fromfile(f, dtype=float64, count=3*N)
    P, chi, B12 = body.reshape((N, 3)).T
    return t, N, P, chi, B12


# function that joins data from all dumps with given number
def rd(no, dump_dir="dumps/"):
    t, N, P, chi, B12 = 0, 0, zeros(0), zeros(0), zeros(0)
    for filename in glob(dump_dir + "dump%05d-??" % no):
        tt, tN, tP, tchi, tB12 = rds(filename)
        if (tt != t) and (t != 0):
            print("Warining: time in %s differs from others" % filename)
        t = tt
        N += tN
        P = append(P, tP)
        chi = append(chi, tchi)
        B12 = append(B12, tB12)
    return t, N, P, chi, B12
