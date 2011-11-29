include make-ddl03910.sys

CFLAGS += -I ./include
LFLAGS += -L ./lib

default: p
.PHONY: clean
clean:
	rm -rf *.o

%.o: %.cpp %.h 
	$(CC) -c $< -o $@ $(CFLAGS)
%.o: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

p: p.o lib/atomistic.o lib/io.o
	$(LD) -o $@ $^ $(LFLAGS)

test-regex: test-regex.o
	$(LD) -o $@ $^ $(LFLAGS)

wfextr: wfextr.o
	$(LD) -o $@ $^ $(LFLAGS)

test-read: test-read.o
	$(LD) -o $@ $^ $(LFLAGS)

test-qi-stack: test-qi-stack.o
	$(LD) -o $@ $^ $(LFLAGS)
test-qi: test-qi.o
	$(LD) -o $@ $^ $(LFLAGS)
