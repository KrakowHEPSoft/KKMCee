INCLUDES = 
METASOURCES = AUTO
lib_LIBRARIES = libDizet.a
libDizet_a_SOURCES = bcqcdl5_14.f bkqcdl5_14.f BornV_StartEW.f dizet6_21.f\
	DZface.f m2tcor5_11.f zf514_aux.f
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
#-----------------------------------------------------------------------------
tableYR.all:	Clean $(COMMAND1) input.all
	(cp  ./inputYR.all ./input.data; time $(COMMAND1))
#-----------------------------------------------------------------------------
table.down:	$(COMMAND1) input.down
	(cp  ./input.down ./input.data; time $(COMMAND1); chmod -w table.down )
table.down-debug:	$(COMMAND1)
	(cp  ./input.down./input.data; xldb -I  ${SRCLIST} $(COMMAND1)  )
#-----------------------------------------------------------------------------
table.up:	$(COMMAND1) input.up
	(cp  ./input.up ./input.data; time $(COMMAND1); chmod -w table.up )
table.up-debug:		$(COMMAND1)
	(cp  ./input.up ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
table.botom:	$(COMMAND1) input.botom
	(cp  ./input.botom ./input.data; time $(COMMAND1); chmod -w table.botom )
botom-debug:		$(COMMAND1)
	(cp  ./input.botom ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
table.mu:	$(COMMAND1) input.mu
	(cp ./input.mu ./input.data; time $(COMMAND1); chmod -w table.mu  )
table.mu-debug:		$(COMMAND1)
	(cp ./input.mu ./input.data; xldb -I  ${SRCLIST} $(COMMAND1)  )
#-----------------------------------------------------------------------------
table.tau:	$(COMMAND1) input.tau
	(cp ./input.tau ./input.data; time $(COMMAND1); chmod -w table.tau  )
table.tau-debug:	$(COMMAND1)
	(cp ./input.tau ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
table.nue:	$(COMMAND1) input.nutau
	(cp ./input.nue ./input.data; time $(COMMAND1); chmod -w table.nue  )
table.nue-debug:	$(COMMAND1)
	(cp ./input.nue ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
table.numu:	$(COMMAND1) input.numu
	(cp ./input.numu ./input.data; time $(COMMAND1); chmod -w table.numu  )
table.numu-debug:	$(COMMAND1)
	(cp ./input.numu ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
table.nutau:	$(COMMAND1) input.nutau
	(cp ./input.nutau ./input.data; time $(COMMAND1); chmod -w table.nutau  )
table.nutau-debug:	$(COMMAND1)
	(cp ./input.nutau ./input.data; xldb -I  ${SRCLIST} $(COMMAND1) )
#-----------------------------------------------------------------------------
tables:	table.down table.up table.botom table.mu table.tau table.nue table.numu table.nutau
##############################################################################
#-----------------------------------------------------------------------------
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
#####################################################################
YR:
	(cp -f ./table.down.345pt.YR    ./table.down.current)
	(cp -f ./table.up.345pt.YR      ./table.up.current)
	(cp -f ./table.botom.345pt.YR   ./table.botom.current)
	(cp -f ./table.mu.345pt.YR      ./table.mu.current)
	(cp -f ./table.tau.345pt.YR     ./table.tau.current)
	(cp -f ./table.nue.345pt.YR     ./table.nue.current)
	(cp -f ./table.numu.345pt.YR    ./table.numu.current)
	(cp -f ./table.nutau.345pt.YR   ./table.nutau.current)
# has to be reproduced!!! old version 340pt
#YRibox:
#	(cp -f ./table.down.YRibox    ./table.down.current)
#	(cp -f ./table.up.YRibox      ./table.up.current)
#	(cp -f ./table.botom.YRibox   ./table.botom.current)
#	(cp -f ./table.mu.YRibox      ./table.mu.current)
#	(cp -f ./table.tau.YRibox     ./table.tau.current)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# this is necessary if grid size changed in bornv.h
Clean:	clean
	(cd ../bornv; make clean)
#####################################################################
#####################################################################
