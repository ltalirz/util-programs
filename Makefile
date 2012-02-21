include make-ipazia-login.sys


# Dependencies
INCDIR     = include
LIBDIR     = lib
BASIC      = $(addprefix $(LIBDIR)/, io.o la.o)
ATOMISTIC  = $(addprefix $(LIBDIR)/atomistic/, fundamental.o)
FORMATS    = $(addprefix $(LIBDIR)/formats/, cube.o xyz.o)
MYLIBDEP   = $(BASIC) $(ATOMISTIC) $(FORMATS)


    
# Programs
STMPROGS   = extrapolate sumbias
UTILPROGS  = cubestride


# Test targets are made like: make test/regex
TEST       = fftw fftw-2 stl blitz inherit karma progress
TESTTARGETS   = $(addprefix test/, $(TEST))
# These targets may depend on my library
TESTLIB    = regex qi qi-stack qi-cptime read write la readcp types p
TESTLIBTARGETS   = $(addprefix test/, $(TESTLIB))
TESTMPI    = mpi
TESTMPITARGETS   = $(addprefix test/, $(TESTLIB))


vpath % include

CFLAGS += -I ./$(INCDIR)
#LFLAGS += -L ./lib  # No archives here

default: p

.PHONY: clean
clean:
	find . -type f -name "*.o" -print | xargs rm

$(LIBDIR)/%.o:: $(LIBDIR)/%.cpp $(INCDIR)/%.h 
	$(CC) -c $< -o $@ $(CFLAGS)

%.o:: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

$(STMPROGS): %: stm/%.o $(MYLIBDEP)
	$(LD) -o bin/$@ $^ $(LFLAGS)
$(UTILPROGS): %: util/%.o $(MYLIBDEP)
	$(LD) -o bin/$@ $^ $(LFLAGS)

# Binaries of test targets are put into test/
$(TESTTARGETS): %: %.o
	$(LD) -o $@ $^ $(LFLAGS)
$(TESTLIBTARGETS): %: %.o $(MYLIBDEP)
	$(LD) -o $@ $^ $(LFLAGS)
