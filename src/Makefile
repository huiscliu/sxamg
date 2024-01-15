
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

CPPFLAGS   := ${CPPFLAGS} ${BUILD_CPPFLAGS} ${USER_CPPFLAGS}
CFLAGS     := ${CFLAGS} ${BUILD_CFLAGS} ${USER_CFLAGS}
SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJS_C  = $(patsubst %.c,%.o,$(DIR_C))

OBJS     := $(OBJS_C)

lib: lib$(SXAMGLIB)$(LIB_SUFFIX)

lib$(SXAMGLIB).a: ${OBJS}
	${AR} -ruv $@ $(OBJS)
	${RANLIB} $@

lib$(SXAMGLIB).so: ${OBJS}
	${BUILD_SHARED} ${OBJS} ${BUILD_SHARED_LIBS} -o libsxamg${LIB_SUFFIX}

.PHONY: clean

clean:
	@rm -fv lib*.a lib*.so
	@rm -fv *.o 
