# Makefile

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O2

# Paths
SRCDIR = src
BUILDDIR = build
TARGET = sav_call

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))

# Include directories
INCLUDES = -I$(SRCDIR)

# Libraries
LIBS = -lhts

# Coverage flags
COVERAGE_FLAGS = -g -O0 --coverage

# Default target
all: $(TARGET)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Build target
$(TARGET): $(BUILDDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

# Build test target
test/$(TARGET): $(BUILDDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

# Compile source files to object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILDDIR)/*.o $(TARGET) build/*.gcda build/*.gcno coverage.info out

# Build with coverage
coverage: CXXFLAGS += $(COVERAGE_FLAGS)
coverage: clean test/$(TARGET)

# Run tests and generate coverage report
test: coverage
	cd test; bash gene_test_data.sh
	cd test; lcov --capture --directory ../build --output-file coverage.info
	cd test; genhtml coverage.info --output-directory out
	@echo "Coverage report generated in 'test/out' directory."
	@echo "Check 'test/check*.err' and 'test/check*.err' to see results on sim data."
