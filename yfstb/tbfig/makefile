COMMAND = tbfig.exe
COMMAND_OBJECTS = tbfig.o
EXTER_LIB2  = ../../glibk/glibk.a
EXTER_SRC2  = ../../glibk/glibk.f 
EXTER_LIB3  = ../../yfslib/yfslib.a
EXTER_SRC3  = ../../yfslib/yfslib.f  
############
all:	$(COMMAND)
FFLAGS  =  -O
LDFLAGS =  -O
F77 = f77
LD  = $(F77)
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB2): $(EXTER_SRC2)
	cd ../../glibk; make
$(EXTER_LIB3): $(EXTER_SRC3)
	cd ../../yfslib; make
####
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)   $(EXTER_LIB2)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB2)
run1:	all
	(./tbfig.exe < ../tbwork1/tbwork1.data )
run2:	all
	(./tbfig.exe < ../tbwork2/tbwork2.data )
run3:	all
	(./tbfig.exe < ../tbwork3/tbwork3.data )
plot2:	run2
	(latex tbwork2.tex; xdvi   tbwork2)
xprint:
	dvips tbfig -o
	xprint -d 53-th tbfig.ps
print:
	dvips tbfig -o
	lpr tbfig.ps
clean:                                                     
	(rm -f *.o *.exe *~ *.dvi *.log *.aux *.ps dump.hst core)
Clean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux *.ps dump.hst core)
	(rm -f *.dvi; rm -f *.aux; rm -f *.log)
	(cd ../../glibk;  make clean)
	(cd ../../yfslib; make clean)
