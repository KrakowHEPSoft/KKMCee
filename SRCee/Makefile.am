
all_includes :=  -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 -I$(top_srcdir)/tauola-photos  $(HEPMC_INCLUDES) -I/opt/PHOTOS-install/include  $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS)  # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking
######CXXFLAGS = -g -O2 -std=c++11 -lstdc++

METASOURCES = AUTO

PHOTOSLIBS =/opt/PHOTOS-install/lib/libPhotospp.so \
			/opt/PHOTOS-install/lib/libPhotosppHepMC3.so
			

#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libKKee.la libKKfm.la
#///////////////////////////////////////////////////////////////
#   Main KKMCee C++ library
#   ROOT Dictionary 1 file definition
dict1C = KK2fDict.cxx
dict1H = KK2fDict.h
#  Header files:
include_HEADERS1 = 	KKee2f.h KKdbase.h KKdizet.h KKborn.h \
					KKevent.h KKlasa.h KKarLud.h KKqed3.h \
					KKarFin.h \
					KKceex.h KKpart.h KKbvir.h \
					TauPair.h HepFace.h
noinst_HEADERS1 = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(all_includes) $(ROOTCINTHDRS1)
	if [ -f "$(dict1H:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dict1H:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
##---------------------------------------------------------------
#  Local library 1
#  Additional flags:
nodist_libKKee_la_SOURCES = $(dict1C)
libKKee_la_SOURCES =	KKee2f.cxx KKdbase.cxx KKdizet.cxx KKborn.cxx \
						KKevent.cxx KKlasa.cxx KKarLud.cxx KKqed3.cxx \
						KKarFin.cxx \
						KKceex.cxx KKpart.cxx KKbvir.cxx \
						TauPair.cxx HepFace.cxx
libKKee_la_LIBADD =  \
	$(PHOTOSLIBS) \
	$(top_builddir)/tauola-photos/libTauolaPhotos.la \
	-lgfortran \
	$(top_builddir)/MCdev/libMCdev.la 
libKKee_la_LDFLAGS = -L$(ROOTLIBDIR) $(LDFLAGS_HEPMC) -L/opt/PHOTOS-install/lib/
#///////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////
#   Library of Foam internal benchmark
#   ROOT Dictionary 2 file definition
dict2C = KKfmDict.cxx
dict2H = KKfmDict.h
#  Header files:
include_HEADERS2 = KKeeFoam.h
noinst_HEADERS2  = LinkDef2.h
#  ROOT Directory generation:
ROOTCINTHDRS2= $(include_HEADERS2) $(noinst_HEADERS2)
$(dict2C) : $(ROOTCINTHDRS2)
	$(ROOTCINT) -f $(dict2C) -c $(all_includes) $(ROOTCINTHDRS2)
	if [ -f "$(dict2H:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dict2H:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
#---------------------------------------------------------------
#  Local library 2
#  Additional flags:
nodist_libKKfm_la_SOURCES = $(dict2C)
libKKfm_la_SOURCES =  KKeeFoam.cxx
libKKfm_la_LIBADD =  \
	$(PHOTOSLIBS)
libKKfm_la_LDFLAGS = -L$(ROOTLIBDIR)
#///////////////////////////////////////////////////////////////
CLEANFILES = $(dict2C) $(dict2H) $(dict1C) $(dict1H)


#///////////////////////////////////////////////////////////////
#  explicit target in ../MCdev
$(top_builddir)/MCdev/libMCdev.la: $(top_builddir)/MCdev/*.h $(top_builddir)/MCdev/*.cxx
	(cd $(top_builddir)/MCdev/; make)
#/////////////////////////////////////////////////////////////

#/////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
#-------------------------------------------------------------------------
bin_PROGRAMS = TestMini1
#/////////////////////////////////////////////////////////
TestMini1_SOURCES = TestMini1.cxx 
TestMini1_LDADD   = \
			$(PHOTOSLIBS) \
			libKKee.la 
##TestMini1_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
TestMini1_LDFLAGS = -L$(ROOTLIBDIR) -L/opt/PHOTOS-install/lib/ -rpath $(ROOTRPATH) -no-install -g

#/////////////////////////////////////////////////////////
test1-start: TestMini1
	(cd ../ make)
	./TestMini1
