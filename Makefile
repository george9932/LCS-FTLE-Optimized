CXX = g++

# change compiler to g++-6 for macOS since the default g++ in
# macOS is a symlink to clang which doesn't support OpenMP
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	CXX = g++-6
endif

CXXFLAGS = -g -Wall -O2 -std=c++17 -fopenmp

INCLUDES = -I./include

SRC = ./fast_computation/discrete_fast_lcs_computation.cpp

BINDIR = ./bin

TARGET = $(BINDIR)/DiscreteFastComputation

OBJ = $(SRC:.cpp=.o)

# Targets
all: prepare_dirs $(TARGET)

# Build target
$(TARGET): $(OBJ) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) -o $(TARGET)

# Object file rule
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRC) -o $(OBJ)

# Directory creation
$(BINDIR):
	mkdir -p $(BINDIR)

# Create directories target
prepare_dirs:
	mkdir -p ./fast_computation/step_flow_maps
	mkdir -p ./fast_computation/results/ftle

# Clean target
clean: prepare_dirs
	rm -f $(TARGET) $(OBJ)

# Phony targets
.PHONY: all clean prepare_dirs
