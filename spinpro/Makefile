#-------------------------------------------------------------------------------
# Test of double bremsstrahlung
#	make Mu-start
# Test of full spin corelations at Born level
#	make Tau-start
#-------------------------------------------------------------------------------
############
all:	demo-start
########### linux flags ###########
#FFLAGS  =  -static -f +B
#LDFLAGS =  
########### collection of HP  flags ###########
### HP options: +se for 100 or more columns, +B for backslash in strings
###             +K  for global SAVE
#FFLAGS  =  -K +se +B
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
#FFLAGS  =  -qfixed=132 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=132 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable
###########
FFLAGS  =  -qfixed=132 -qextname -O
LDFLAGS =  -qextname -bloadmap:lmap
##
F77 = f77
LD  = $(F77)
#-----------------------------------------------------------------------
SRCLIST= "../KK2f ../glibk ../jetset ../bornv "
#-----------------------------------------------------------------------
###########
EXT_LIB1  = ../KK2f/KK2f.a
EXT_SRC1  = ../KK2f/*.f		../KK2f/*.h
EXT_LIB2  = ../glibk/glibk.a
EXT_SRC2  = ../glibk/*.f	../glibk/*.h
EXT_LIB3  = ../jetset/jetset.a
EXT_SRC3  = ../jetset/*.f
EXT_LIB4  = ../bornv/bornv.a
EXT_SRC4  = ../bornv/*.f	../bornv/*.h
EXT_LIB5  = ../tauola/glib.a
EXT_SRC5  = ../tauola/*.f
EXT_LIB6  = ../photos/glib.a
EXT_SRC6  = ../photos/photos.f
# double bremsstrahlung
EXT_LIB7  = ../2brems/lib/mainlib.a
EXT_SRC7  = ../2brems/lib/*.f
EXT_LIB8  = ../2brems/trzynue/lib/Z0_lib.a
EXT_SRC8  = ../2brems/trzynue/lib/*.f
EXT_LIB9  = ../2brems/trzynue/analysis/analysis_lib.a
EXT_SRC9  = ../2brems/trzynue/analysis/*.f
EXT_LIB10 = ../2brems/2brems.a
EXT_SRC10 = ../2brems/*.f
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXT_LIB1): $(EXT_SRC1)
	(cd ../KK2f; make)
$(EXT_LIB2): $(EXT_SRC2)
	(cd ../glibk; make)
$(EXT_LIB3): $(EXT_SRC3)
	(cd ../jetset; make)
$(EXT_LIB4): $(EXT_SRC4)
	(cd ../bornv; make)
$(EXT_LIB5): $(EXT_SRC5)
	(cd ../tauola; make)
$(EXT_LIB6): $(EXT_SRC6)
	(cd ../photos; make)
$(EXT_LIB7): $(EXT_SRC7)
	(cd ../2brems/lib; make)
$(EXT_LIB8): $(EXT_SRC8)
	(cd ../2brems/trzynue/lib; make)
$(EXT_LIB9): $(EXT_SRC9)
	(cd ../2brems/trzynue/analysis; make)
$(EXT_LIB10): $(EXT_SRC10)
	(cd ../2brems; make)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Link procedure
COMMAND1 = MainPro.exe
OBJECTS1 = MainPro.o
$(COMMAND1): $(OBJECTS1)  \
	$(EXT_LIB1)  $(EXT_LIB2)  $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) \
                 $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6)
#-----------------------------------------------------------------------
# Link procedure
COMMAND2 = ProdMC.exe
OBJECTS2 = ProdMC.o
$(COMMAND2): $(OBJECTS2)  \
	$(EXT_LIB1)  $(EXT_LIB2)  $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7) $(EXT_LIB8) $(EXT_LIB9) $(EXT_LIB10)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS2) \
                 $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7) $(EXT_LIB8) $(EXT_LIB9) $(EXT_LIB10)
#-----------------------------------------------------------------------
Tau-start:	$(COMMAND1)
	(cd ./Tau;  make start )
Tau-stop:
	(cd ./Tau;  make stop )
Tau-start-debug:	$(COMMAND1)
	(cd ./Tau;  make start-debug )
Tau-cont-debug:	$(COMMAND1)
	(cp Tau/pro.hst.110 Tau/pro.hst )
	(cp Tau/semaphore.110 Tau/semaphore )
###	(cd ./Tau;  make cont )
	(cd ./Tau;  make cont-debug )
#-----------------------------------------------------------------------
Mu-start:	$(COMMAND2)
	(cd ./Mu;  make start )
Mu-stop:
	(cd ./Mu;  make stop )
Mu-cont:
	(cd ./Mu;  make cont )
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
#-----------------------------------------------------
	(sed -e 's#<main>#MainPro#' -e 's#<dset>#Tau#' -e 's#<class>#long#'\
		./makefile.production > ./Tau/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Mu#' -e 's#<class>#long#'\
		./makefile.production > ./Mu/makefile)
#-----------------------------------------------------
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
#=============================================================================
# Update fortran compilation flags
makflag:
	(cd ../KK2f;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../jetset;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../tauola;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../photos;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../glibk;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../dizet;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../bornv;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../bornv/tabtest;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../2brems;			sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../2brems/lib;		sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../2brems/trzynue/lib;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../2brems/trzynue/analysis;	sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
#-----------------------------------------------------------------------
Clean:
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd ../glibk;	make clean)
	(cd ../KK2f;	make clean)
	(cd ../jetset;	make clean)
	(cd ../tauola;	make clean)
	(cd ../photos;	make clean)
	(cd ../bornv;	make clean)
	(cd ../dizet;	make clean)
	(cd ../bornv/tabtest;	make clean)
	(cd ../2brems;			make clean)
	(cd ../2brems/lib;		make clean)
	(cd ../2brems/trzynue/lib;	make clean)
	(cd ../2brems/trzynue/analysis;	make clean)
#-----------------------------------------------------------------------

