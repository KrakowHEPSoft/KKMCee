COMMAND = tbmain.exe
COMMAND_OBJECTS = tbmain.o erichter.o syost.o
EXTER_LIB1  = ../yfs3/yfs3.a
EXTER_SRC1  = ../yfs3/yfs3mod.f ../yfs3/yfs3.f
EXTER_LIB2  = ../glibk/glibk.a
EXTER_SRC2  = ../glibk/glibk.f  
EXTER_LIB3  = ../yfslib/yfslib.a
EXTER_SRC3  = ../yfslib/yfslib.f 
############
all:	$(COMMAND)
FFLAGS  =  -O 
LDFLAGS =  -O 
F77 = f77
LD  = $(F77)
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB1): $(EXTER_SRC1)
	cd ../yfs3; make
$(EXTER_LIB2): $(EXTER_SRC2)
	cd ../glibk; make
$(EXTER_LIB3): $(EXTER_SRC3)
	cd ../yfslib; make
#############
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB1)  $(EXTER_LIB2) $(EXTER_LIB3)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) \
              $(EXTER_LIB1) $(EXTER_LIB2) $(EXTER_LIB3)
run1:	all
	(cd ./tbwork1; ./go)
run2:	all
	(cd ./tbwork2; ./go)
run3:	all
	(cd ./tbwork3; ./go)
plot2:	
	(cd ./tbfig; make plot2)
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd tbfig; make clean)
Clean:
	(rm -f *.o *.a *~ *.exe *.dvi *.log *.aux core)
	(cd tbfig;     make Clean)
	(cd ../yfslib; make clean)
	(cd ../glibk;  make clean)
	(cd ../yfs3;   make clean)
fortex:
	tex syost.f; xdvi syost

