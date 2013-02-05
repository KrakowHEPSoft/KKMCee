COMMAND = phodem.exe
COMMAND_OBJECTS = phodem.o
EXTER_LIB1  = ../../jetset/jetset.a
EXTER_SRC1  = ../../jetset/jetset74.f 
EXTER_LIB2  = ../glib.a
EXTER_SRC2  = ../photos.f
############
all:	$(COMMAND)
#FFLAGS  = -static
#LDFLAGS = 
########### collection of HP  flags ###########
### HP options: +se for 100 or more columns, +B for backslash in strings
###             +K  for global SAVE
#FFLAGS  =  -K +se +B
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
#FFLAGS  =  -qfixed=100 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=100 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable
###########
FFLAGS  =  -qfixed=100 -qextname -O
LDFLAGS =  -qextname -bloadmap:lmap

F77 = f77
LD  = $(F77)
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB1): $(EXTER_SRC1)
	cd ../../jetset; make
$(EXTER_LIB2): $(EXTER_SRC2)
	cd ../; make
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB1) $(EXTER_LIB2)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB1) $(EXTER_LIB2)
#$(COMMAND): $(COMMAND_OBJECTS)  $(EXTER_LIB1) 
#	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS) $(EXTER_LIB1) 
run:
	cd prod; sh go
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe; rm -f *.log; rm -f *.ps; rm -f *.dvi; rm -f *.aux
