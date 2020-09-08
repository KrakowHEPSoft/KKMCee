
METASOURCES = AUTO
lib_LIBRARIES = libDizet.a
libDizet_a_SOURCES = DZface.f BornV_StartEW.f  hadr5n95.f hadr5n16.f hadr5n17.f zfmai6_42.f zfusr6_44.f zfbib6_40.f dizet6_45.f APV_lib.f acol6_1p.f m2tcor5_11.f pairho6_40.f funang6_30.f expifi6_30.f bcqcdl5_14.f bkqcdl5_14.f bhang4_640.f zf630_aux.f 
bin_PROGRAMS = TabMain
TabMain_SOURCES = TabMain.f
TabMain_LDADD = $(top_builddir)/dizet/libDizet.a


##########################################################################
##############################################################################
COMMAND1 = ./TabMain
SRCLIST= "./"
#-----------------------------------------------------------------------------
table.all:	Clean $(COMMAND1) input.all
	(cp  ./input.all ./input.data; time $(COMMAND1))
#	make table-clean	make table.all		make table-check
table-clean:
	(rm -f  table.tau.current table.down.current table.up.current table.botom.current table.mu.current \
		table.nue.current table.numu.current table.nutau.current)
table-check:
	(diff table.down.current   table.down    | more)
	(diff table.up.current     table.up      | more)
	(diff table.botom.current  table.botom   | more)
	(diff table.mu.current     table.mu      | more)
	(diff table.tau.current    table.tau     | more)
	(diff table.nue.current    table.nue     | more)
	(diff table.numu.current   table.numu    | more)
	(diff table.nutau.current  table.nutau   | more)
#####################################################################
default:
	(cp -f ./table.down.345pt.default    ./table.down.current)
	(cp -f ./table.up.345pt.default      ./table.up.current)
	(cp -f ./table.botom.345pt.default   ./table.botom.current)
	(cp -f ./table.mu.345pt.default      ./table.mu.current)
	(cp -f ./table.tau.345pt.default     ./table.tau.current)
	(cp -f ./table.nue.345pt.default     ./table.nue.current)
	(cp -f ./table.numu.345pt.default    ./table.numu.current)
	(cp -f ./table.nutau.345pt.default   ./table.nutau.current)
#####################################################################
table-cp:
	(cp table.down   table.down.current    )
	(cp table.up     table.up.current      )
	(cp table.botom  table.botom.current   )
	(cp table.mu     table.mu.current      )
	(cp table.tau    table.tau.current     )
	(cp table.nue    table.nue.current     )
	(cp table.numu   table.numu.current    )
	(cp table.nutau  table.nutau.current   )
table-cp-default:
	(cp table.down   table.down.345pt.default    )
	(cp table.up     table.up.345pt.default      )
	(cp table.botom  table.botom.345pt.default   )
	(cp table.mu     table.mu.345pt.default      )
	(cp table.tau    table.tau.345pt.default     )
	(cp table.nue    table.nue.345pt.default     )
	(cp table.numu   table.numu.345pt.default    )
	(cp table.nutau  table.nutau.345pt.default   )
	(chmod -w table.*.345pt.default)
table-cp-YR:
	(cp table.down   table.down.345pt.YR    )
	(cp table.up     table.up.345pt.YR      )
	(cp table.botom  table.botom.345pt.YR   )
	(cp table.mu     table.mu.345pt.YR      )
	(cp table.tau    table.tau.345pt.YR     )
	(cp table.nue    table.nue.345pt.YR     )
	(cp table.numu   table.numu.345pt.YR    )
	(cp table.nutau  table.nutau.345pt.YR   )
table-cp-YRibox:
	(cp table.down   table.down.YRibox    )
	(cp table.up     table.up.YRibox      )
	(cp table.botom  table.botom.YRibo    )
	(cp table.mu     table.mu.YRibox      )
	(cp table.tau    table.tau.YRibox     )
#-----------------------------------------------------------------------------
# this is necessary if grid size changed in bornv.h
Clean:	clean
	(cd ../bornv; make clean)
#####################################################################
#####################################################################
AM_FFLAGS = @MY_F77_FLAGS@
