if HAVE_GMP_WITH_CXX
SUBDIRS = gmp
endif

EXTRA_DIST =																	\
	333 333.gro.chk 333.lat 333.mar														\
	333.1 333.1.gro.chk 333.1.mar 333.1.weights 333.1.weights.max										\
	333.2 333.2.gro.chk 333.2.mar 333.2.weights 333.2.weights.max										\
	333.3 333.3.gro.chk 333.3.mar 333.3.weights 333.3.weights.max										\
	334 334.gro.chk 334.lat 334.mar														\
	335 335.gro.chk 335.lat 335.mar														\
	344 344.gro.chk 344.lat 344.mar														\
	cuww1 cuww1.1 cuww1.1.mar cuww1.1.sign cuww1.cost cuww1.gro.chk cuww1.lat cuww1.mar							\
	cuww2 cuww2.cost cuww2.gro.chk cuww2.lat cuww2.mar											\
	cuww2.trunc1 cuww2.trunc1.cost cuww2.trunc1.gro.chk cuww2.trunc1.lat cuww2.trunc1.mar cuww2.trunc1.weights cuww2.trunc1.weights.max	\
	cuww2.trunc2 cuww2.trunc2.cost cuww2.trunc2.gro.chk cuww2.trunc2.lat cuww2.trunc2.mar cuww2.trunc2.weights cuww2.trunc2.weights.max	\
	cuww3 cuww3.cost cuww3.gro.chk cuww3.lat cuww3.mar											\
	cuww4 cuww4.cost cuww4.gro.chk cuww4.lat cuww4.mar											\
	cuww5 cuww5.cost cuww5.gro.chk cuww5.lat cuww5.mar											\
	grin grin.cost grin.gro.chk grin.lat grin.mar												\
	grin10900.1 grin10900.1.cost grin10900.1.gro.chk grin10900.1.mar grin10900.1.zsol							\
	grin10900.2 grin10900.2.cost grin10900.2.gro.chk grin10900.2.mar grin10900.2.zsol							\
	grin10900.3 grin10900.3.cost grin10900.3.gro.chk grin10900.3.mar grin10900.3.zsol							\
	grin10900.4 grin10900.4.cost grin10900.4.gro.chk grin10900.4.mar grin10900.4.zsol							\
	grin10900.5 grin10900.5.cost grin10900.5.mar grin10900.5.zsol


check-fifo: 
	@../check.template $(srcdir) "-fifo" src/groebner/groebner -q -afifo

check-weighted:
	@../check.template $(srcdir) "-weighted" src/groebner/groebner -q -aweighted

check-local: check-fifo check-weighted
	@../flavorsum "-fifo" "-weighted"

CLEANFILES = .num_failures .num_failures-fifo .num_failures-weighted


