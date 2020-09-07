
#F77     = g77
F77 = gfortran
FFLAGS  = -fno-automatic -g -ffixed-line-length-132 -fno-backslash
CC      = g++
CFLAGS  = -fPIC -Wall -Wno-write-strings
LD      = $(F77)
LDFLAGS = -fno-automatic -g -lstdc++
AR      = ar
RANLIB  = ranlib


CPP     = /lib/cpp 
CPPFLAGS= -C -P -traditional
DIFF    = diff -ib
