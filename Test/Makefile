.POSIX:
.SUFFIXES:
CXX      = clang++
#LCGR    = LCG_96/x86_64-slc6-gcc8-opt# Now define this in bashrc
SHELL    = /bin/sh

TARGET   = dptest
OBJD     = build
SRCD     = src
SRCS     = $(shell /bin/ls  $(SRCD)/*.cpp)
OBJS     = $(shell for n in $(SRCS) ; do echo $(OBJD)/`basename $$n .cpp`.o ; done)
DEPS     = $(shell for n in $(SRCS) ; do echo $(OBJD)/`basename $$n .cpp`.d ; done)

#INCD   != find $(SRCD) -type d
#INCF   !=$(addprefix -I,$(INCD))# -isystem${CPATH}
INCF     =-I$(SRCD) -I/cvmfs/sft.cern.ch/lcg/views/$(LCGR)/include
#CXFLAGS         = $(shell if test "clang++" = $(CXX) ; then echo "-Weverything"\
#       "-Wno-c++98-compat -Wno-double-promotion -Wno-covered-switch-default" ; fi)
CXFLAGS +=$(INCF) -MMD -MP -std=c++17 -march=native -pipe -O2 -g\
        -Wall -Wextra -Wpedantic -mfma# -flto -fPIE
#ROOTDIR         = $(shell root-config --libdir)
ROOTLIB  = $(shell root-config --libs  )
#LDFLAGS         =-L$(ROOTDIR)
LDFLAGS +=$(ROOTLIB) -L/cvmfs/sft.cern.ch/lcg/views/$(LCGR)/lib\
        -Wl,--enable-new-dtags#,-pie -fuse-ld=gold

$(OBJD)/$(TARGET) :             $(OBJS)
	$(CXX) $(CXFLAGS) -o $(OBJD)/calchisto.o        -c TestDeep.cxx
	$(CXX) $(CXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
$(OBJD)/%.o : $(SRCD)/%.cpp
	@mkdir -p $(OBJD)
	$(CXX) $(CXFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	-rm    $(OBJD)/$(TARGET) $(OBJD)/*.o $(OBJD)/*.d
	-rmdir $(OBJD)

-include $(DEPS)

