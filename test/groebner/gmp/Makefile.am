EXTRA_DIST = 																					\
	335.gro.chk 335.mar 335.mat																		\
	cuww1 cuww1.cols cuww1.cost cuww1.gro.chk cuww1.lat cuww1.mar														\
	cuww1.1 cuww1.1.mar cuww1.1.sign cuww2.trunc1 cuww2.trunc1.cost cuww2.trunc1.gro.chk cuww2.trunc1.lat cuww2.trunc1.mar cuww2.trunc1.weights cuww2.trunc1.weights.max	\
	grin10900.1 grin10900.1.cost grin10900.1.gro.chk grin10900.1.mar grin10900.1.zsol


check-fifo: 
	@../../check.template $(srcdir) "-fifo" src/groebner/groebner -q -p arb -afifo

check-weighted:
	@../../check.template $(srcdir) "-weighted" src/groebner/groebner -q -p arb -aweighted

check-local: check-fifo check-weighted
	@../../flavorsum "-fifo" "-weighted"

CLEANFILES = .num_failures .num_failures-fifo .num_failures-weighted

