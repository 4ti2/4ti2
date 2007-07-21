EXTRA_DIST = 												\
	344 344.cost 344.cost.start 344.gro.chk 344.gro.start 344.mar					\
	cuww1 cuww1.cost cuww1.cost.start cuww1.gro.chk cuww1.gro.start					\
	grin grin.cost grin.cost.start grin.gro.chk grin.gro.start grin.mar				\
	prob20 prob20.cost prob20.cost.start prob20.gro.chk prob20.gro.start prob20.lat prob20.sign

check-local:
	@../check.template $(srcdir) "" src/groebner/walk -q 

CLEANFILES = .num_failures
