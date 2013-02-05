COMMAND = trzynue_test.exe
COMMAND_OBJECTS = trzynue_test.o
EXTER_LIB1  = lib/Z0_lib.a
EXTER_SRC1  = lib/Z0_test2ph.f lib/brelib.f lib/dublib.f lib/declib.f lib/trzylib.f lib/Z0_test3ph.f
EXTER_LIB2  = ../lib/mainlib.a
EXTER_SRC2  = ../lib/glibnew.f ../lib/photos_20.f
EXTER_LIB3  = analysis/analysis_lib.a
EXTER_SRC3  = analysis/bukery.f 
############
all:	$(COMMAND)
FFLAGS  = -static -O
LDFLAGS = -static -O
F77 = f77
LD  = $(F77)
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXTER_LIB1): $(EXTER_SRC1)
	cd  lib; make
$(EXTER_LIB2): $(EXTER_SRC2)
	cd ../lib; make
$(EXTER_LIB3): $(EXTER_SRC3)
	cd analysis; make
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)   $(EXTER_LIB1) $(EXTER_LIB2) $(EXTER_LIB3)
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS)  $(EXTER_LIB1) $(EXTER_LIB2) $(EXTER_LIB3)
run:
	cd prod; sh go
clean:                                                     
	rm -f *.o *.a  *~ *.exe
Clean:	clean
	cd ../lib; rm -f *.o *.a  *~ *.exe
	cd lib;    rm -f *.o *.a  *~ *.exe
	cd analysis;    rm -f *.o *.a  *~ *.exe
	cd ploty;    rm -f *.o *.a  *~ *.exe
gtar:   Clean
	(cd ../..; tar cvf amplitudy.tar amplitudy)
	(cd ../..; gzip amplitudy.tar)
uu: gtar
	(cd ../..; uuencode amplitudy.tar.gz amplitudy.tar.gz > amplitudy.tar.gz.uu)


