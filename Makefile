include make-ddl03910.sys

default: p
.PHONY: clean
clean:
	rm -rf *.o

%.o: %.cpp %.h 
	$(CC) -o $@ -c $(CFLAGS) $<
%.o: %.cpp
	$(CC) -o $@ -c $(CFLAGS) $<

p.x: p.o include/classes.o
	$(LD) -o $@ $^ $(LFLAGS)

test-parse.x: test-parse.o
	$(LD) -o $@ $^ $(LFLAGS)

wfextr.x: wfextr.o
	$(LD) -o $@ $^ $(LFLAGS)
