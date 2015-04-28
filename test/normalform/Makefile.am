EXTRA_DIST = \
	grin grin.cost grin.feas grin.gro grin.nf.chk \
	4coins.mat 4coins.cost 4coins.feas 4coins.gro 4coins.nf.chk 

check-local:
	@../check.template $(srcdir) "" src/groebner/normalform -q 

CLEANFILES = .num_failures
