include make-ddl03910.sys

CFLAGS += -I ./include
LFLAGS += -L ./lib

vpath % include:lib

default: p
.PHONY: clean
clean:
	find . -type f -name "*.o" -print | xargs rm

%.o: %.cpp %.h 
	$(CC) -c $< -o $@ $(CFLAGS)
%.o: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

p: p.o atomistic.o io.o
	$(LD) -o $@ $^ $(LFLAGS)

test-regex: test-regex.o
	$(LD) -o $@ $^ $(LFLAGS)

wfextr: wfextr.o atomistic.o io.o
	$(LD) -o $@ $^ $(LFLAGS)

test-read: test-read.o
	$(LD) -o $@ $^ $(LFLAGS)

test-qi-stack: test-qi-stack.o
	$(LD) -o $@ $^ $(LFLAGS)

test-qi: test-qi.o
	$(LD) -o $@ $^ $(LFLAGS)
test-blitz: test-blitz.o
	$(LD) -o $@ $^ $(LFLAGS)
