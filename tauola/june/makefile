COMMAND = taumain.exe
COMMAND_OBJECTS = taumain.o
EXTER_LIB1  = ../glib.a
EXTER_SRC1  = ../formf.f ../tauola.f 
EXTER_LIB2  = ../../glibk/glibk.a
EXTER_SRC2  = ../../glibk/GLK.f 
EXTER_LIB3  = ../../jetset/jetset.a
EXTER_SRC3  = ../../jetset/jetset74.f 
EXTER_LIB4  = ../../photos/glib.a
EXTER_SRC4  = ../../photos/photos.f
############
all:	$(COMMAND)
#FFLAGS  = -K -O 
#LDFLAGS = -K -O 
########### collection of HP  flags ###########
### HP options: +se for 100 or more columns, +B for backslash in strings
###             +K  for global SAVE
#FFLAGS  =  -K +se +B
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
FFLAGS  =  -qfixed=100 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=100 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable
###########
#FFLAGS  =  -qfixed=100 -qextname -O
LDFLAGS =  -qextname -bloadmap:lmap

F77 = f77
LD  = $(F77)
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB1): $(EXTER_SRC1)
	cd ..; make
$(EXTER_LIB2): $(EXTER_SRC2)
	cd ../../glibk; make
$(EXTER_LIB3): $(EXTER_SRC3)
	cd ../../jetset; make
$(EXTER_LIB4): $(EXTER_SRC4)
	cd ../../photos; make
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB1)  $(EXTER_LIB2) $(EXTER_LIB3)  $(EXTER_LIB4)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB1) $(EXTER_LIB2) $(EXTER_LIB3)  $(EXTER_LIB4)
#$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB1) 
#	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB1) 
run:
	cd prod1; sh go
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
