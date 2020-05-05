.POSIX:
.SUFFIXES:
CXX   ?= g++
SHELL  = /bin/sh

TARGET = eta
OBJD   = build
SRCD   = src
SRCS   = $(wildcard $(SRCD)/*.cpp)
OBJS   = $(patsubst $(SRCD)/%.cpp,$(OBJD)/%.cpp.o,$(SRCS))
DEPS   = $(OBJS:.o=.d)

INCD    =$(shell find $(SRCD) -type d)
INCFLAGS=$(addprefix -I,$(INCD)) -isystem${CPATH}

CXXFLAGS=$(INCFLAGS) -MMD -MP -std=c++17 -march=native -pipe\
	-O3 -Wall -Wextra -Wpedantic -Wold-style-cast

LDPATH =-L$(shell root-config --libdir) \
	-L/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/lib
LIBS   =  $(shell root-config --libs)
LDFLAGS=-rdynamic ${LDPATH} ${LIBS} \
	    -Wl,-Rlib,-R../lib,-R${PWD}/lib,--enable-new-dtags

$(OBJD)/$(TARGET) :        $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

$(OBJD)/ :
	mkdir -p $(OBJD)

$(OBJD)/%.cpp.o : $(SRCD)/%.cpp | $(OBJD)/
	$(CXX) $(CXXFLAGS)  -c $< -o $@

.PHONY: clean

clean:
	-rm    $(OBJD)/$(TARGET) $(OBJD)/*.o $(OBJD)/*.d
	-rmdir $(OBJD)

-include $(DEPS)
