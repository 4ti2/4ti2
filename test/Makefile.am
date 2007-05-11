SUBDIRS =					\
	normalform				\
	qsolve					\
	graver					\
	gensymm					\
	groebner				\
	zsolve					\
	walk					\
	ppi					\
	zbasis					\
	output					\
	markov					\
	circuits				\
	genmodel				\
	minimize				\
	rays					\
	hilbert

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
