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
.PRECIOUS: all essential extended special

#entry targets
all: all.res $(GECCOPRG)
	rmdir $(SCRDIR)

essential: essential.res $(GECCOPRG)
	rmdir $(SCRDIR)

extended:  extended.res $(GECCOPRG)
	rmdir $(SCRDIR)

special:   special.res $(GECCOPRG)
	rmdir $(SCRDIR)


#targets to run the test (defined in Makefile.parallel) , print summary and clean up

export SCRDIR 
export GECCOPRG

all.res: $(GECCOPRG)
	@rm -f *.err ;\
	$(MAKE) -f Makefile.parallel all;\
	echo "Result for whole suite:" ;\
	. ./message.sh ;\
	print_all *.err |tee all.res ;\
	print_summary *.err | tee -a all.res ;\

essential.res: $(GECCOPRG)
	@rm -f *.e.err ;\
	$(MAKE) -f Makefile.parallel essential;\
	echo "Result for essential test set:" ;\
	. ./message.sh ;\
	print_all *.e.err |tee essential.res ;\
	print_summary *.e.err | tee -a essential.res ;\


extended.res:   $(GECCOPRG)
	@rm -f *.x.err;\
	$(MAKE) -f Makefile.parallel extended;\
	echo "Result for extended test set:" ;\
	. ./message.sh ;\
	print_all *.x.err |tee extended.res ;\
	print_summary *.x.err | tee -a extended.res ;\


special.res:  $(GECCOPRG) 
	@rm -f *.s.err ;\
	$(MAKE) -f Makefile.parallel special;\
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

