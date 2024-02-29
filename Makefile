# Copyright (C) Panua Technologies, 2023

# Set C++ compiler and compiler flags
CXX = g++
CXX_FLAGS = -O3

ARCH:=$(shell uname -m)
ifeq ($(ARCH),x86_64)
    # x86 specific flags
	PARDISO_INSTALL_DIR=$(HOME)/local/panua-pardiso-x86
	TARGET_SUFFIX=x86
else ifeq ($(ARCH),aarch64)
    # ARM specific flags
	PARDISO_INSTALL_DIR=$(HOME)/local/panua-pardiso-arm
	TARGET_SUFFIX=arm
endif

# Additional libraries
LIBS=\
-lm\
-L$(PARDISO_INSTALL_DIR)/lib -lpardiso -Wl,--rpath=$(PARDISO_INSTALL_DIR)/lib\

INCFLAGS=\
-I$(PARDISO_INSTALL_DIR)/include\
-I$(HOME)/files/fmca\
-I$(HOME)/local/eigen/include/eigen3\
-I pardiso_interface.h\

# Source directory
SRC_DIR = .

# Source files - Automatically detect all .cpp and .f90 files in the source directory
CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
CPP_OBJS = $(CPP_SRCS:.cpp=.o)

# Targets
all: test

test: $(CPP_OBJS)
	@echo ">> Pardiso Location: $(PARDISO_INSTALL_DIR)"
	@echo ">> Linking $@"
	$(CXX) $(CXX_FLAGS) $(INCFLAGS) $^ $(LIBS) -o $@.exe_$(TARGET_SUFFIX)

%.o: %.cpp
	@echo ">> Compiling $@"
	$(CXX) $(CXX_FLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -f *.o

purge:
	rm -f *.o *.exe*