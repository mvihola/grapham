include ../Makefile.in
GRAPHAM_ROOT=$(shell cd ..; pwd)
COPTS:=-O3 -Wall -std=c99 -pedantic -fshort-enums -fPIC \
-I$(GRAPHAM_ROOT)/src
LDOPTS=-L$(GRAPHAM_ROOT)/src -lgrapham -lm

OBJS:=dbanana.so dcircular.so dcross.so thresholds.so

default: $(OBJS)

%.o: %.c
	$(CC) $(COPTS) $(CFLAGS) -c -o $@ $<

%.so: %.o
	$(CC) -shared -o $@ $< $(LDOPTS)

.PHONY: clean
clean:
	rm -f $(OBJS) 
