# Compiler to be used
CC = g++

# Compiler flags
CFLAGS = -Wall -std=c++17 -pthread -Ofast -s
# INCLUDES = -I ~/mambaforge/include 
LIBS = -lyaml-cpp
# LIBS = -L ~/mambaforge/lib -Wl,--disable-new-dtags,-rpath,~/mambaforge/lib -llapack -lblas -larmadillo -lyaml-cpp

# Source and build directories
SRCDIR = ./src
BUILDDIR = ./obj
TARGETDIR = ./bin

# Application name (target)
TARGET = LifNet

# Get all source files
SRCS := $(wildcard $(SRCDIR)/*.cpp)

# Object files are the same as sources but with .o extension
OBJS := $(SRCS:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Build both Release and Debug by default
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGETDIR)/$@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)/*.o $(TARGETDIR)/$(TARGET)

# Recompile if any headers change
.PHONY: clean
