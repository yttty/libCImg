# Terrill Yang
# Mar 10, 2017

CXX ?= clang++
# CFLAGS = -L/usr/X11R6/lib
CFLAGS = -I../../CImg -Iinclude -g
LFLAGS = -lm -lpthread -lX11

CC_SOURCES = $(shell find . -name "*.cc")
OBJECTS = $(patsubst %.cc, %.o, $(CC_SOURCES))
BINARY = edge

$(BINARY) : $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $(BINARY)

$(OBJECTS) : %.o: %.cc
	$(CXX) -c $(CFLAGS) $< -o $@

.PHONY : clean 
clean : 
	rm $(BINARY) $(OBJECTS)
