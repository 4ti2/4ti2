EXTRA_DIST =					\
	333.mat 333.mar.chk			\
	cuww1.mat cuww1.mar.chk			\
	cuww3.mat cuww3.mar.chk			\
	cuww5.mat cuww5.mar.chk

check-local:
	@../../check.template $(srcdir) "" src/groebner/markov -q -gmax-min

CLEANFILES = .num_failures
