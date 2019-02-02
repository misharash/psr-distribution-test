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

# function that reads all data from dump
from numpy import fromfile, float64


def rd(no):
    f = open("dumps/dump%05d" % no, "rb")
    # read header
    header = f.readline().decode()
    L = header.split()
    t = float(L[0])
    N = int(L[1])
    # read body
    body = fromfile(f, dtype=float64, count=3*N)
    P, chi, B12 = body.reshape((N, 3)).T
    return t, N, P, chi, B12
