
METASOURCES = AUTO
lib_LIBRARIES = libDizet.a libDizet2.a
libDizet_a_SOURCES = \
	DZface.f BornV_StartEW.f \
	dizet6_42.f zf630_aux.f bcqcdl5_14.f bkqcdl5_14.f m2tcor5_11.f hadr5n17_compact.f 

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
	dizet6_42.f zf630_aux.f bcqcdl5_14.f bkqcdl5_14.f m2tcor5_11.f hadr5n17_compact.f 
