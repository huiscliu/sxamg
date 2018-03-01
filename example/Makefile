
include ../Makefile.inc

default: lib  amg amg-low-level

lib:
	@(cd ../; make)

amg.o: amg.c ../src/lib$(NLIB).a
amg-low-level.o: amg-low-level.c ../src/lib$(NLIB).a

clean:
	rm -f amg amg.o amg-low-level.o amg-low-level
