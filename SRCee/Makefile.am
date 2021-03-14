
all_includes :=  -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking
######CXXFLAGS = -g -O2 -std=c++11 -lstdc++

METASOURCES = AUTO

#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libKKee.la
#///////////////////////////////////////////////////////////////
#   Main KKMChh C++ library
#   ROOT Dictionary 1 file definition
dict1C = KK2fDict.cxx
dict1H = KK2fDict.h
CLEANFILES = $(dict1C) $(dict1H)
#  Header files:
include_HEADERS1 = 	KKee2f.h KKdbase.h KKdizet.h KKborn.h \
					KKevent.h KKlasa.h KKarLud.h KKqed3.h
#					KKarFin.h \
#                   KKceex.h KKpart.h KKbvir.h
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
						KKevent.cxx KKlasa.cxx KKarLud.cxx KKqed3.cxx
#						KKarFin.cxx \
#						KKceex.cxx KKpart.cxx KKbvir.cxx
libKKee_la_LIBADD =  $(top_builddir)/MCdev/libMCdev.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libKKee_la_LDFLAGS = -L$(ROOTLIBDIR)
#///////////////////////////////////////////////////////////////


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
			$(top_builddir)/SRCee/libKKee.la \
			$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
TestMini1_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g

#/////////////////////////////////////////////////////////
test1-start: TestMini1
	(cd ../ make)
	./TestMini1
