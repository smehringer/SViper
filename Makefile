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
# or download into current folder and set to:
ifndef SEQAN_LIB
$(error [ATTENTION] Variable SEQAN_LIB is not set. Please export the inlcude path)
endif

CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1

# set sviper inlcude path
CXXFLAGS+=-I$(INCLUDEDIR)

# linked libraries
LDLIBS=-lz -lpthread

# enable openMP
CXXFLAGS+=-fopenmp

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
DEBUGCXXFLAGS+=$(CXXFLAGS)
DEBUGCXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
RELEASECXXFLAGS+=$(CXXFLAGS)
RELEASECXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DNDEBUG

all: sviper

debug: $(SRCDIR)/sviper.cpp
	 $(CXX) $(DEBUGCXXFLAGS) -o sviper_debug $(SRCDIR)/sviper.cpp $(LDLIBS)

sviper: $(SRCDIR)/sviper.cpp
	 $(CXX) $(RELEASECXXFLAGS) -o sviper $(SRCDIR)/sviper.cpp $(LDLIBS)

utilities: $(TARGETS)

$(UTILITYDIR)/%: $(SRCDIR)/%.cpp
	@mkdir -p $(UTILITYDIR)
	$(CXX) $(RELEASECXXFLAGS) $(LDLIBS) -o $@ $<

clean:
	rm -f sviper
	rm -f sviper_debug
	@echo Cleaning executables in utilities
	$(shell find utilities -type f -regex '^[^.]+' -delete)

.PHONY: clean

