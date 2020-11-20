GRAPHAM_ROOT=$(shell cd ..; pwd)
include ../Makefile.in

ifeq ($(GRAPHAM_VERSION),)
	GRAPHAM_VERSION:=$(shell date +'%Y%m%d')
endif

CFLAGS+=-D_GRAPHAM_VERSION=$(GRAPHAM_VERSION) -D_GRAPHAM_ROOT=$(GRAPHAM_ROOT)
# Dunno what -fomit-frame-pointer and -fshort-enums do, but they
# messed up living with luamatrix! 
COPTS+=-O3 -Wall -std=c99 -pedantic -fshort-enums -fPIC
# This will clear a warning message, but will produce slower code.
#COPTS:=-O3 -Wall -std=c99 -pedantic -fshort-enums -fno-strict-referecing
LIBS:=-Wl,-E -ldl -lm

LIB:=libgrapham.a
OBJS:=grapham.o grapham_math.o lua_tools.o grapham_io.o grapham_core.o \
custom_type.o vector_type.o grapham_lib.o grapham_rand.o lua_math.o 
LIB_OBJS:=grapham_math.o
LDOPTS+=$(GRAPHAM_ROOT)/netlib/netlib.a

ifeq ($(ENVIRONMENT), mingw)
	EXECNAME=grapham.exe
else
	EXECNAME=grapham
endif

ifneq ($(DSFMT_PATH),)
	CFLAGS+=-D_HAVE_DSFMT -I$(DSFMT_PATH)
	OBJS+=$(DSFMT_PATH)/dSFMT.o
endif

ifneq ($(LUA_PATH),)
	CFLAGS+=-I$(LUA_PATH)/include -I$(LUA_PATH)/include/lua5.1 \
-I$(LUA_PATH)/src
	LDOPTS+=-L$(LUA_PATH)/lib -I$(LUA_PATH)/lib/lua/5.1 \
-L$(LUA_PATH)/src
	LIBS+=-llua
else
	CFLAGS+=-I/usr/include/lua5.1 -I/usr/local/include/lua5.1
	LDOPTS+=-L/usr/lib/lua/5.1 -L/usr/local/lib/lua/5.1
	LIBS+=-llua5.1
endif

ifneq ($(NUMLUA_PATH),)
	CFLAGS+=-D_HAVE_NUMLUA -I$(NUMLUA_PATH)/src
	LDOPTS+=$(NUMLUA_PATH)/src/luamatrix.a $(NUMLUA_PATH)/src/luacomplex.a \
$(NUMLUA_PATH)/src/luaspfun.a $(NUMLUA_PATH)/src/luarng.a \
-L$(NUMLUA_PATH)/lib -lfnlib -ldcdflib  -llapack -lblas -lranlib -lg2c
	OBJS+=matrix_type.o
endif

default: $(EXECNAME) lib
lib: $(LIB) 
$(LIB): $(LIB_OBJS)
	cp $(GRAPHAM_ROOT)/netlib/netlib.a $@
	ar ru $@ $(LIB_OBJS)

$(DSFMT_PATH)/dSFMT.o: $(DSFMT_PATH)/dSFMT.c
	cd "$(DSFMT_PATH)"; make dSFMT.o

%.o: %.c 
	$(CC) $(COPTS) $(CFLAGS) -c -o $@ $<

$(EXECNAME): $(OBJS) 
	$(CC) -o $@ $(OBJS) $(LDOPTS) $(LIBS) 

.PHONY: clean lib default
clean:
	rm -f $(EXECNAME) $(OBJS) $(LIB)
