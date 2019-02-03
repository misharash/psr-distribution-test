#!/bin/bash
all:
	g++ -Wall -O3 main.cpp create.cpp approx.cpp dump.cpp evol.cpp delete.cpp birth.cpp restart.cpp
