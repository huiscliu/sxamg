
default: lib  amg amg-low-level

lib:
	@(cd ../; make)

include ../Makefile.inc

amg.o: amg.c ../lib/lib$(SXAMGLIB).a
amg-low-level.o: amg-low-level.c ../lib/lib$(SXAMGLIB).a

clean:
	rm -f amg amg.o amg-low-level.o amg-low-level
