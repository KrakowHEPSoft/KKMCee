AC_INIT([KKMCee],[4.25])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([configure.ac])
# ':='-style assignments are not portable
AM_INIT_AUTOMAKE([-Wall -Wno-portability])
AC_PREFIX_DEFAULT(${PWD}/)
LT_INIT
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_CC

AC_PROG_F77([gfortran f95])

AM_CONFIG_HEADER(config.h)

# Do compilation tests using the compiler, preprocessor, and file extensions
AC_LANG([Fortran 77])
###AC_LANG_FORTRAN77

AC_F77_LIBRARY_LDFLAGS

MY_F77_FLAGS="-O -g -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash -fPIC"
AC_SUBST(MY_F77_FLAGS)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Uncomment only for making distro to be used without automake
###AM_MAINTAINER_MODE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ROOT_PATH(3.02/00, , AC_MSG_ERROR( ROOT has not been found.))
#HAS_CERNLIB( , AC_MSG_ERROR( CERNLIB has not been found))
AC_LANG(C)
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([[#include<stdio.h>]],
			[[const char hw[] = "Hello, World\n";]],
                       [[fputs (hw, stdout);]])],
	[AC_MSG_RESULT(TESTPROG OK)], [AC_MSG_RESULT(TestProg Failed)]
	);

##### macros defining different queue and class 
##### for NQS batch scripts
QUEUE_CMD
CLASS_CMD 


AC_CONFIG_FILES([Makefile \
dizet/Makefile \
jetset/Makefile \
tauola-photos/Makefile  \
Foam2/Makefile \
MCdev/Makefile MCdev/farming/Makefile \
SRCee/Makefile \
DZface/Makefile \
ProdRun/Makefile  ProdRun/workMini/Makefile  ProdRun/work1/Makefile ProdRun/workFoam/Makefile \
                  ProdRun/workNU/Makefile ProdDigest/Makefile \
])
AC_OUTPUT