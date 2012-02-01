include make-ipazia-login.sys

CFLAGS += -I ./include
#LFLAGS += -L ./lib  # No archives here

MYLIBDIR   = lib
MYLIBOBJ   = io.o la.o atomistic.o
MYLIBDEP   = $(addprefix $(MYLIBDIR)/, $(MYLIBOBJ))

STMPROGS   =  extrapolate

# Test targets are made like: make test/regex
TESTP  =  p regex qi-stack qi karma
TESTP +=  fftw fftw-2 stl blitz
TESTP +=  read write la
TESTTARGETS   = $(addprefix test/, $(TESTP))

vpath % include

default: p

.PHONY: clean
clean:
	find . -type f -name "*.o" -print | xargs rm

$(MYLIBDIR)/%.o:: %.cpp %.h 
	$(CC) -c $< -o $@ $(CFLAGS)

%.o:: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

$(STMPROGS): %: stm/%.o $(MYLIBDEP)
	$(LD) -o bin/$@ $^ $(LFLAGS)

# Binaries of test targets are put into test/
$(TESTTARGETS): %: %.o $(MYLIBDEP)
	$(LD) -o $@ $^ $(LFLAGS)
