# Terrill Yang
# Apr 5, 2017

CXX = clang++
CFLAGS = -I../.. -g -Wall -std=c++11
CFLAGS += -Dcimg_use_openmp -O3
# CFLAGS += -L/usr/X11R6/lib
LFLAGS = -lm -lpthread -lX11
LFLAGS += -fopenmp

CC_SOURCES = $(shell find . -name "*.cc")
OBJECTS = $(patsubst %.cc, %.o, $(CC_SOURCES))
LIBCIMG = ../../libcimg.o
BINARY = color_transfer

$(BINARY) : $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBCIMG) $(LFLAGS) -o $(BINARY)

$(OBJECTS) : %.o: %.cc
	make -C ../..
	$(CXX) -c $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm $(BINARY) $(OBJECTS)
