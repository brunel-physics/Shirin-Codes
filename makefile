CXX = g++

LDPATH := -L$(shell root-config --libdir) \
		 -L/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/lib
LIBS := $(shell root-config --libs)
LDFLAGS := -rdynamic ${LDPATH} ${LIBS} \
		  -Wl,-Rlib,-R../lib,-R${PWD}/lib,--enable-new-dtags

TARGET_EXEC ?= eta

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS)) -isystem${CPATH}

CXXFLAGS ?= $(INC_FLAGS) -MMD -MP -std=c++17 -march=native -mtune=native -pipe \
			-O3 -Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual \
			-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
			-Wlogical-op -Wmissing-declarations -Wmissing-include-dirs \
			-Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls \
			-Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel \
			-Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused \
			-Wuseless-cast -Wzero-as-null-pointer-constant -Wduplicated-cond \
			-Wduplicated-branches -Wrestrict -Wnull-dereference \
			-Wdouble-promotion -Wimplicit-fallthrough=4

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p
