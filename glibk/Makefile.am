
METASOURCES = AUTO

lib_LIBRARIES = libGlibk.a
libGlibk_a_SOURCES = GLK.f MathLib.f PltLog.f PseuMar.f ReaData.f
bin_PROGRAMS = combine
combine_SOURCES = combine.f
combine_LDADD = $(top_builddir)/glibk/libGlibk.a

AM_FFLAGS = @MY_F77_FLAGS@
