
default: lib  amg amg-low-level cg-amg cg-amg2 gmres-amg

lib:
	@(cd ../; make)

include ../Makefile.inc

amg.o: amg.c ../src/lib$(SXAMGLIB)$(LIB_SUFFIX)
amg-low-level.o: amg-low-level.c ../src/lib$(SXAMGLIB)$(LIB_SUFFIX)
cg-amg.o: cg-amg.c ../src/lib$(SXAMGLIB)$(LIB_SUFFIX)
cg-amg2.o: cg-amg2.c ../src/lib$(SXAMGLIB)$(LIB_SUFFIX)
gmres-amg.o: gmres-amg.c ../src/lib$(SXAMGLIB)$(LIB_SUFFIX)

clean:
	rm -f *.o amg amg-low-level cg-amg cg-amg2 gmres-amg
