# By Terrill Yang
# libCimg - Utilities for the Cimg library

CXX ?= clang++
CFLAGS = -ICImg -g -Wall -std=c++11
CFLAGS += -Dcimg_use_openmp -O3
# CFLAGS += -L/usr/X11R6/lib

CC_SOURCES = libcimg.cc
OBJECTS = $(patsubst %.cc, %.o, $(CC_SOURCES))

$(OBJECTS) : %.o: %.cc
	$(CXX) -c $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm $(OBJECTS)
