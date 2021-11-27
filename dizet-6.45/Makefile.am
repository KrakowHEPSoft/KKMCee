
METASOURCES = AUTO
lib_LIBRARIES = libDizet.a libDizet2.a

libDizet_a_SOURCES = \
	DZface.f BornV_StartEW.f  \
	hadr5n95.f hadr5n16.f hadr5n17.f zfmai6_42.f zfusr6_44.f zfbib6_40.f \
	dizet6_45.f APV_lib.f acol6_1p.f m2tcor5_11.f pairho6_40.f funang6_30.f expifi6_30.f \
	bcqcdl5_14.f bkqcdl5_14.f bhang4_640.f zf630_aux.f 

bin_PROGRAMS = TabMain
TabMain_SOURCES = TabMain.f
TabMain_LDADD = $(top_builddir)/dizet/libDizet.a
#####################################################################
COMMAND1 = ./TabMain
#--------------------------------------------------------------------
EWtables:	$(COMMAND1) input.all
	(cp  ./input.all ./input.data; time $(COMMAND1))
#####################################################################
AM_FFLAGS = @MY_F77_FLAGS@

libDizet2_a_SOURCES = \
	hadr5n95.f hadr5n16.f hadr5n17.f zfmai6_42.f zfusr6_44.f zfbib6_40.f \
	dizet6_45.f APV_lib.f acol6_1p.f m2tcor5_11.f pairho6_40.f funang6_30.f expifi6_30.f \
	bcqcdl5_14.f bkqcdl5_14.f bhang4_640.f zf630_aux.f
