SRCDIR := src
BUILDDIR := build
UTILITYDIR := utilities
INCLUDEDIR := include

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name utilities*.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TARGETS := $(patsubst $(SRCDIR)/%,$(UTILITYDIR)/%,$(SOURCES:.$(SRCEXT)=))

# Use version > 4.9.2 of g++
CXX=g++
CC=$(CXX)

CXXFLAGS+=-std=c++14

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB=/nfs/prog/bioinfo/apps-x86_64/seqan-library/2.2.0/include

CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1
CXXFLAGS+=-I./src/
CXXFLAGS+=-I$(INCLUDEDIR)
LDLIBS=-lz -lpthread
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
#CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DNDEBUG

all: sviper evaluate_final_mapping

sviper: $(SRCDIR)/sviper.cpp
	 $(CXX) $(CXXFLAGS) $(LDLIBS) -fopenmp -o sviper $(SRCDIR)/sviper.cpp

evaluate_final_mapping: $(SRCDIR)/evaluate_final_mapping.cpp
	 $(CXX) $(CXXFLAGS) $(LDLIBS) -fopenmp -o evaluate_final_mapping $(SRCDIR)/evaluate_final_mapping.cpp

utilities: $(TARGETS)

$(UTILITYDIR)/%: $(SRCDIR)/%.cpp
	@mkdir -p $(UTILITYDIR)
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $@ $<

clean:
	rm sviper
	rm evaluate_final_mapping
	@echo Cleaning executables in utilities
	$(shell find utilities -type f -regex '^[^.]+' -delete)

.PHONY: clean

