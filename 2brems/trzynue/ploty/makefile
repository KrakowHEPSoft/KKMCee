COMMAND = dallas.exe
COMMAND_OBJECTS = dallas.o
EXTER_LIB3  = ../../lib/mainlib.a
EXTER_SRC3  = ../../lib/glibnew.f ../../lib/photos_20.f
#######################################
all:	$(COMMAND)
#FFLAGS  = -K -O
#LDFLAGS = -K -O
FFLAGS  = -O
LDFLAGS = -O
F77 = f77
LD = $(F77)
#######################################
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB3): $(EXTER_SRC3)
	cd ../../lib; make
$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB3)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB3)
run:
	(./dallas.exe)
clean:
	(rm -f *.o; rm -f *.exe; rm -f *~; rm -f *.log; rm -f *.aux; rm -f *.ps)
	(rm -f dump.hst)
view:
	(cd output; latex dallas; xdvi dallas)
	(cd ..)
print:
	(cd output)
	(dvips dallas.dvi -o  dallas.ps)
	(ghostview dallas.ps)
	(cd ..)