include make.sel

##### Flags
INCDIR     = include
LIBDIR     = lib

vpath % $(INCDIR)

CFLAGS += -I ./$(INCDIR)
#LFLAGS += -L ./lib  # No archives here

##### Dependencies
BASIC      = io la
ATOMISTIC  = $(addprefix atomistic/, fundamental)
FORMATS    = $(addprefix formats/, cp cp2k cube espresso xyz gnuplot stm)
WRAPPERS   = $(addprefix wrappers/, lapack)

COMPONENTS = $(BASIC) $(ATOMISTIC) $(FORMATS) $(WRAPPERS)
INCDEP     = $(addprefix $(INCDIR)/, $(addsuffix .hpp, $(COMPONENTS)))
LIBDEP     = $(addprefix $(LIBDIR)/, $(addsuffix .o, $(COMPONENTS)))

##### Programs
STMPROGS   = extrapolate sumbias sts stm extrapolate2
UTILPROGS  = cubestride cubescale cubesquare cuberoot cubeabs cubezprofile
UTILPROGS += cubediravg espressowfn cubeplane cubeprod


# Test targets are made like: make test/regex
TEST       = fftw fftw-2 stl blitz inherit karma progress po core lapack
TESTTARGETS   = $(addprefix test/, $(TEST))
# These targets may depend on my library
TESTLIB    = regex qi qi-stack qi-cptime read write la readcp types p stm readesp
TESTLIB   += lap wrap
TESTLIBTARGETS   = $(addprefix test/, $(TESTLIB))
TESTMPI    = mpi
TESTMPITARGETS   = $(addprefix test/, $(TESTLIB))

##### Targets

default: p

.PHONY: clean
clean:
	find . -type f -name "*.o" -print | xargs rm

$(LIBDIR)/%.o:: $(LIBDIR)/%.cpp $(INCDIR)/%.hpp 
	$(CC) -c $< -o $@ $(CFLAGS)
# General programs might depend on headers of the library
%.o:: %.cpp $(INCDEP)
	$(CC) -c $< -o $@ $(CFLAGS)

$(STMPROGS): %: stm/%.o $(LIBDEP)
	$(LD) -o bin/$@ $^ $(LFLAGS)
$(UTILPROGS): %: util/%.o $(LIBDEP)
	$(LD) -o bin/$@ $^ $(LFLAGS)

# Binaries of test targets are put into test/
$(TESTTARGETS): %: %.o
	$(LD) -o $@ $^ $(LFLAGS)
$(TESTLIBTARGETS): %: %.o $(LIBDEP)
	$(LD) -o $@ $^ $(LFLAGS)
