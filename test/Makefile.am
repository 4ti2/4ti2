SUBDIRS =					\
	gensymm					\
	ppi					\
	output					\
	genmodel

if ENABLE_GROEBNER
SUBDIRS +=						\
	circuits					\
	groebner					\
	markov						\
	minimize					\
	rays						\
	walk						\
	qsolve						\
	zbasis						\
	normalform				
endif

if ENABLE_ZSOLVE			
SUBDIRS +=					\
	graver					\
	zsolve					\
	hilbert
endif

# The check target runs automatic checks to verify that 4ti2 was compiled
# correctly and is running correctly.

# Override the check-recursive directory, so all subdirs are check in parallel
# (if -j is used).  This can save some real time.
check-recursive: $(SUBDIRS:%=checkdir-%)
	$(MAKE) check-local 

checkdir-%: %
	cd $< && $(MAKE) check

check-local:
	@./subdirsum $(SUBDIRS)
	@./reportfailures

CLEANFILES = .num_failures
