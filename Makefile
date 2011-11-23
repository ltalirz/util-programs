include make-ipazia-login.sys

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

test-parse: test-parse.o
	$(LD) -o $@.x $^ $(LFLAGS)

wfextr: wfextr.o
	$(LD) -o $@.x $^ $(LFLAGS)
