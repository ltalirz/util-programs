CC = g++
LD = g++
CFLAGS = -I ~/apps/boost/boost_1_48_0/
#CFLAGSROB = -I /opt/intel/ict_2011/mkl/include/fftw -I ~/apps/roberto-progs/ -Wno-deprecated
LFLAGS = -L ~/apps/boost/boost_1_48_0/binaries/lib/ -lboost_program_options -static
#LFLAGSROB = -L ~/apps/roberto-progs/ -L /opt/intel/ict_2011/mkl/lib/intel64/  -lfftw3xf_intel -lmkl_intel_lp64   -lmkl_intel_thread -lmkl_core -L /opt/intel/ict_2011/composerxe-2011.0.084/compiler/lib/intel64 -liomp5 -lpthread -lm

APPDIR = /home/tal127/apps

default: p
.PHONY: clean
clean:
	rm -rf *.o

%.o: %.cpp %.h 
	$(CC) -c $(CFLAGS) $<
%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

p: p.o
	$(LD) -o $@.x $^ $(LFLAGS)

wfextr: wfextr.o
	$(LD) -o $@.x $^ $(LFLAGS)
