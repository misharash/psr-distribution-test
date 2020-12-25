#!/bin/bash
CXXFLAGS = -Wall -O3 -MMD
SRC = $(wildcard *.cpp)
OBJS = ${SRC:.cpp=.o}
DEPS = ${OBJS:.o=.d}
EXEC = psr_dist
.PHONY: all clean
all: $(EXEC)
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
clean:
	rm -f *.o *.d $(EXEC)
-include ${DEPS}
