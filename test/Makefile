#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#



topdir:=..
srcdir:=.

include $(topdir)/make.config

SCRDIR:=${GECCO_TMP}

#include $(topdir)/make.rules

GECCOPRG:=$(CURDIR)/$(topdir)/gecco/$(bindir)/gecco.x

.NOTPARALLEL:

.PHONY: all essential special extended res_clean out_clean err_clean clean essential.res special.res extended.res all.res
.PRECIOUS: all essential extended special delscr 




#entry targets

%.clean :
	@ rm -f $(basename $@).err $(basename $@).out $(basename $@).res

%.test: $(GECCOPRG) | %.clean
	@-rm -f $(basename $@).err
	@$(MAKE) -f Makefile --no-print-directory mkscr
	@-$(MAKE) -f Makefile.parallel --no-print-directory $(basename $@).res
	@-$(MAKE) -f Makefile --no-print-directory delscr
	@. ./message.sh ;\
	print_all $(basename $@).err ;\
	print_summary $(basename $@).err 

all:  $(GECCOPRG)
	@$(MAKE) -f Makefile --no-print-directory mkscr
	@-$(MAKE) -f Makefile --no-print-directory all.res
	@-$(MAKE) -f Makefile --no-print-directory delscr

essential:  $(GECCOPRG)
	@$(MAKE) -f Makefile --no-print-directory mkscr
	@-$(MAKE) -f Makefile --no-print-directory essential.res
	@-$(MAKE) -f Makefile --no-print-directory delscr

extended:  $(GECCOPRG)
	@$(MAKE) -f Makefile --no-print-directory mkscr
	@-$(MAKE) -f Makefile --no-print-directory extended.res
	@-$(MAKE) -f Makefile --no-print-directory delscr

special:   $(GECCOPRG)
	@$(MAKE) -f Makefile --no-print-directory mkscr
	@-$(MAKE) -f Makefile --no-print-directory special.res
	@-$(MAKE) -f Makefile --no-print-directory delscr


#targets to run the test (defined in Makefile.parallel) , print summary and clean up


essential_tests:=$(wildcard *.e.chk)
extended_tests:=$(wildcard *.x.chk)
special_tests:=$(wildcard *.s.chk)
essential_outputs:=$(addsuffix .out, $(basename $(essential_tests)))
extended_outputs:=$(addsuffix .out, $(basename $(extended_tests)))
special_outputs:=$(addsuffix .out, $(basename $(special_tests)))
essential_results:=$(addsuffix .res, $(basename $(essential_tests)))
extended_results:=$(addsuffix .res, $(basename $(extended_tests)))
special_results:=$(addsuffix .res, $(basename $(special_tests)))


export SCRDIR 
export GECCOPRG

all.res: $(GECCOPRG)
	@rm -f *.err ;\
	$(MAKE) -f Makefile.parallel $(essential_results) $(extended_results) $(special_results);\
	echo "Result for whole suite:" ;\
	. ./message.sh ;\
	print_all *.err |tee all.res ;\
	print_summary *.err | tee -a all.res ;\

essential.res: $(GECCOPRG)
	@rm -f *.e.err ;\
	$(MAKE) -f Makefile.parallel $(essential_results);\
	echo "Result for essential test set:" ;\
	. ./message.sh ;\
	print_all *.e.err |tee essential.res ;\
	print_summary *.e.err | tee -a essential.res ;\


extended.res:   $(GECCOPRG)
	@rm -f *.x.err;\
	$(MAKE) -f Makefile.parallel $(extended_results);\
	echo "Result for extended test set:" ;\
	. ./message.sh ;\
	print_all *.x.err |tee extended.res ;\
	print_summary *.x.err | tee -a extended.res ;\


special.res:  $(GECCOPRG) 
	@rm -f *.s.err ;\
	$(MAKE) -f Makefile.parallel $(special_results);\
	echo "Result for special test set:" ;\
	. ./message.sh ;\
	print_all *.s.err |tee special.res ;\
	print_summary *.s.err | tee -a special.res ;\


#cleaning targets

err_clean:
	@set -e ; rm -f *.err

res_clean:
	rm -f *.res ;  rm -f *.sum

out_clean:
	rm -f *.out 

clean:	err_clean res_clean out_clean

delscr:
	@rmdir $(SCRDIR) || echo "could not remove scratch directory: "$(SCRDIR)" please remove by hand."

mkscr:       # is removed in the toplevel targets in the other Makefile
	@if [ ! -d $(SCRDIR) ];\
	then\
	  mkdir $(SCRDIR);\
	fi;

