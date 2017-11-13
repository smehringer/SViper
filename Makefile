SRCDIR := src
BUILDDIR := build
TARGETDIR := utilities

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TARGETS := $(patsubst $(SRCDIR)/%,$(TARGETDIR)/%,$(SOURCES:.$(SRCEXT)=))

# Use version > 4.9.2 of g++
CXX=g++
CC=$(CXX)

CXXFLAGS+=-std=c++14

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB=/nfs/prog/bioinfo/apps-x86_64/seqan-library/2.2.0/include

CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1
CXXFLAGS+=-I./src/
LDLIBS=-lz -lpthread -fopenmp
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
#CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DNDEBUG

all: $(TARGETS)

$(TARGETDIR)/%: $(SRCDIR)/%.cpp
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $@ $<

clean:
	@echo Cleaning executables in utilities
	$(shell find utilities -type f -regex '^[^.]+' -delete)

.PHONY: clean

