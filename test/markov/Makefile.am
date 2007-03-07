SUBDIRS = saturation
if HAVE_GMP_WITH_CXX
SUBDIRS += gmp
endif

EXTRA_DIST = 									\
	333 333.mar.chk								\
	334 334.mar.chk								\
	335 335.lat 335.mar.chk							\
	cuww1 cuww1.bs cuww1.mar.chk cuww1.vs					\
	cuww2 cuww2.mar.chk							\
	cuww3 cuww3.A cuww3.mar.chk						\
	cuww4 cuww4.mar.chk							\
	cuww5 cuww5.mar.chk							\
	grin grin.mar.chk							\
	grin10900.1.cost grin10900.1.lat grin10900.1.mar.chk grin10900.1.zsol	\
	grin10900.2.cost grin10900.2.lat grin10900.2.mar.chk grin10900.2.zsol	\
	grin10900.3.cost grin10900.3.lat grin10900.3.mar.chk grin10900.3.zsol	\
	grin10900.4.cost grin10900.4.lat grin10900.4.mar.chk grin10900.4.zsol	\
	grin1412 grin1412.mar.chk						\
	prob02 prob02.mar.chk							\
	prob04 prob04.mar.chk							\
	prob06 prob06.lat prob06.mar.chk					\
	prob08 prob08.mar.chk							\
	prob10 prob10.mar.chk							\
	prob12 prob12.mar.chk							\
	prob14 prob14.mar.chk							\
	prob16 prob16.lat prob16.mar.chk					\
	prob18 prob18.mar.chk							\
	prob20 prob20.mar.chk

check-hybrid:
	@../check.template $(srcdir) "-hybrid" src/groebner/markov -q -ghybrid

check-project-and-lift:
	@../check.template $(srcdir) "-project-and-lift" src/groebner/markov -q -gproject-and-lift

check-local: check-hybrid check-project-and-lift
	@../flavorsum "-hybrid" "-project-and-lift"
	@../subdirsum $(SUBDIRS)

CLEANFILES = .num_failures .num_failures-hybrid .num_failures-project-and-lift
