AUTOMAKE_OPTIONS = foreign 1.6
SUBDIRS = glibk bornv dizet KK2f jetset tauola photos KKsem ffbench MaMar
###SUBDIRS = glibk bornv dizet KK2f jetset tauola photos KKsem yfspro MaMar ffbench

# Automake 1.9.1 manual suggtests putting all local-macrros in m4 subdir
# and using this:
ACLOCAL_AMFLAGS = -I m4

#Trick: we want to distribute file called "makefile" in /work1.
#We want to force the use of "Makefile" (instead of "makefile", which
#has a precedence) in the whole installation process.
AM_MAKEFLAGS = -f Makefile

DATE = 2012-06
VERS = v5.00
