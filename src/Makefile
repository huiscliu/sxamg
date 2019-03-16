
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJS_C  = $(patsubst %.c,%.o,$(DIR_C))

OBJS     := $(OBJS_C)

lib: lib$(SXAMGLIB).a

lib$(SXAMGLIB).a: ${OBJS}
	${AR} -ruv $@ $(OBJS)
	${RANLIB} $@

.PHONY: clean

clean:
	@rm -fv lib*.a
	@rm -fv *.o 
