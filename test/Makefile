#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

topdir:=..
srcdir:=.

include $(topdir)/make.config

# some new rules (must be defined here BEFORE the objects)
%.res : %.out
	@set -e ; $(CHECK) $(CHECKFLAGS) $(basename $<)

%.out : %.chk $(GECCOPRG)
	@set -e ; $(RUNTEST) $(RUNFLAGS) $(basename $<)

SCRDIR=


GECCOPRG:=$(CURDIR)/$(topdir)/gecco/$(bindir)/gecco.x

RUNTEST:=./runtest.py
CHECK:=./runtest.py
ifndef SCRDIR
  RUNFLAGS:=-r -p $(GECCOPRG)
else
  RUNFLAGS:=-r -p $(GECCOPRG) -s $(SCRDIR)
endif
CHECKFLAGS:=-c

include $(topdir)/make.rules

essential_tests:=$(wildcard *.e.chk)
extended_tests:=$(wildcard *.x.chk)
special_tests:=$(wildcard *.s.chk)
essential_outputs:=$(addsuffix .out, $(basename $(essential_tests)))
extended_outputs:=$(addsuffix .out, $(basename $(extended_tests)))
special_outputs:=$(addsuffix .out, $(basename $(special_tests)))
essential_results:=$(addsuffix .res, $(basename $(essential_tests)))
extended_results:=$(addsuffix .res, $(basename $(extended_tests)))
special_results:=$(addsuffix .res, $(basename $(special_tests)))


.PHONY: all essential special extended res_clean out_clean err_clean clean message

all: essential extended special $(GECCOPRG)

essential: essential.res $(GECCOPRG)

extended:  extended.res $(GECCOPRG)

special:   special.res $(GECCOPRG)

essential.res: $(essential_results) $(essential_outputs) $(GECCOPRG)
	@set -e ; echo "Result for essential test set:" ;\
	$(MAKE) -s message | tee essential.res

extended.res:  $(extended_results)  $(extended_outputs)  $(GECCOPRG)
	@set -e ; echo "Result for extended test set:" ;\
	$(MAKE) -s message | tee extended.res

special.res:   $(special_results)   $(special_outputs) 
	@set -e ; echo "Result for special test set:" ;\
	$(MAKE) -s message | tee special.res

err_files:=$(wildcard *buggy) $(wildcard *failed)

message:
	@set -e ; \
	if [ -f check_OK ] ; then \
          echo "============================================="; \
	  echo " Successful tests:" ; \
          echo "---------------------------------------------"; \
	  cat check_OK ; \
          echo "============================================="; \
	  rm -f check_OK ; \
	fi
	@set -e ; \
	if [ -f check_failed ] ; then \
          echo "============================================="; \
	  echo "Failed tests:" ; \
          echo "---------------------------------------------"; \
	  cat check_failed ; \
          echo "============================================="; \
	  rm -f check_failed ; \
	fi
	@set -e ; \
	if [ -f check_buggy ] ; then \
          echo "============================================="; \
	  echo "The following tests have buggy checks:" ; \
          echo "---------------------------------------------"; \
	  cat check_buggy ; \
          echo "============================================="; \
	  rm -f check_buggy ; \
	fi
	@set -e ; \
	if [ -f run_failed ] ; then \
          echo "============================================="; \
	  echo "The following tests did not run properly:" ; \
          echo "---------------------------------------------"; \
	  cat run_failed ; \
          echo "============================================="; \
	  rm -f run_failed ; \
	fi
	@set -e ; \
	if [ -f setup_buggy ] ; then \
          echo "============================================="; \
	  echo "The following tests had buggy inputs:" ; \
          echo "---------------------------------------------"; \
	  cat setup_buggy ; \
          echo "============================================="; \
	  rm -f setup_buggy ; \
	fi
	@set -e ; \
	if [ -z "$(err_files)" -o "$(err_files)" = " " ] ; then \
	  echo "A lucky day: all tests seem to run correctly!" ; \
	else  \
	  echo "Oh no: some tests give errors, see above!" ; \
	  echo "!!! ERRORS !!!" ; \
        fi

err_clean:
	@set -e ; rm -f check_OK check_failed check_buggy run_OK run_failed setup_buggy

res_clean:
	rm -f *.res 

out_clean:
	rm -f *.out 

clean:	err_clean res_clean out_clean

