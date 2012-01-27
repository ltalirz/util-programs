include make-ipazia-login.sys

CFLAGS += -I ./include
#LFLAGS += -L ./lib  # No archives here

MYLIBDIR   = lib
MYLIBOBJ   = io.o la.o atomistic.o
MYLIBDEP   = $(addprefix $(MYLIBDIR)/, $(MYLIBOBJ))

MYPROGS    = p wfextr test-read test-write test-la

# These programs should not depend on my library
TESTPROGS  =  test-regex test-qi-stack test-qi test-karma
TESTPROGS +=  test-fftw test-fftw-2 test-stl test-blitz
TESTPROGS +=  

vpath % include

default: p

.PHONY: clean
clean:
	find . -type f -name "*.o" -print | xargs rm

$(MYLIBDIR)/%.o:: %.cpp %.h 
	$(CC) -c $< -o $@ $(CFLAGS)

%.o:: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

$(TESTPROGS): %: %.o
	$(LD) -o $@ $^ $(LFLAGS)

$(MYPROGS): %: %.o $(MYLIBDEP)
	$(LD) -o $@ $^ $(LFLAGS)


