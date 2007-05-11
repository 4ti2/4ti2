EXTRA_DIST = 									\
	335 335.mar.chk								\
	cuww3 cuww3.A cuww3.mar.chk						\
	grin10900.3.cost grin10900.3.lat grin10900.3.mar.chk grin10900.3.zsol	\
	grin1412 grin1412.mar.chk						\
	prob06 prob06.lat prob06.mar.chk

check-hybrid:
	@../../check.template $(srcdir) "-hybrid" src/groebner/markov -q -p arb -ghybrid

check-project-and-lift:
	@../../check.template $(srcdir) "-project-and-lift" src/groebner/markov -q -p arb -gproject-and-lift

check-local: check-hybrid check-project-and-lift
	@../../flavorsum "-hybrid" "-project-and-lift"

CLEANFILES = .num_failures .num_failures-hybrid .num_failures-project-and-lift
