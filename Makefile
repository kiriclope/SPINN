# Compiler to be used
CC = g++

# Compiler flags
CFLAGS = -Wall -std=c++17 -Ofast -s -march=native -funroll-loops -ftree-vectorize -ffast-math -fomit-frame-pointer -fexpensive-optimizations
CFLAGS_DEBUG = -Wall -Wextra -std=c++17 -O0 -g

# INCLUDES = -I ~/home/leon/models/lif_cpp/include
# LIBS = -I/home/leon/homebrew/Cellar/yaml-cpp/0.8.0/include -L/home/leon/homebrew/Cellar/yaml-cpp/0.8.0/lib -lyaml-cpp
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

debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGETDIR)/$@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)

# Create a special rule for the debug build
$(TARGET)_debug: $(OBJS)
	$(CC) -o $(TARGETDIR)/$@ $^ $(CFLAGS_DEBUG) $(INCLUDES) $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
ifndef DEBUG
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
else
	$(CC) $(CFLAGS_DEBUG) $(INCLUDES) $(LIBS) -c $< -o $@
endif

clean:
	rm -rf $(BUILDDIR)/*.o $(TARGETDIR)/$(TARGET)

# Recompile if any headers change
.PHONY: clean debug
