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

# HTSlib paths (adjust if necessary)
HTSLIB_INC = /usr/local/include
HTSLIB_LIB = /usr/local/lib

# Adjust these paths if HTSlib is installed elsewhere
CXXFLAGS += -I$(HTSLIB_INC)
LDFLAGS = -L$(HTSLIB_LIB)

# Default target
all: $(TARGET)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Build target
$(TARGET): $(BUILDDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

# Compile source files to object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILDDIR)/*.o $(TARGET)
