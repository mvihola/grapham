include ../Makefile.in

FOPTS+=-O3 -Wall -fPIC

LIB:=netlib.a
LIB_OBJS:=daxpy.o dchdc.o dchdd.o dchud.o ddot.o dlgama.o dnrm2.o drotg.o \
dscal.o dswap.o dtrdi.o

default: lib
lib: $(LIB) 
$(LIB): $(LIB_OBJS)
	ar rcs $@ $(LIB_OBJS)

%.o: %.f
	$(FC) $(FOPTS) -c -o $@ $<

.PHONY: clean lib default
clean:
	rm -f $(OBJS) $(LIB) $(LIB_OBJS)
