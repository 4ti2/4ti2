SUBDIRS =
if HAVE_GMP_WITH_CXX
SUBDIRS += gmp
endif

EXTRA_DIST = 								\
	33.1 33.1.qhom.chk 33.1.sign					\
	33.2 33.2.qhom.chk 33.2.sign					\
	33.3 33.3.qhom.chk 33.3.sign					\
	33.c 33.c.qhom.chk 33.c.sign					\
	33.r 33.r.qhom.chk 33.r.sign					\
	44.c.A 44.c.lat 44.c.qhom.chk 44.c.sign				\
	44.r.A 44.r.lat 44.r.qhom.chk 44.r.sign				\
	55.r 55.r.qhom.chk 55.r.sign					\
	ex ex.qhom.chk ex.rel ex.sign					\
	ppp.r ppp.r.qfree.chk ppp.r.qhom.chk ppp.r.sign			\
	small small.qhom.chk small.rel small.sign			\
	subspace subspace.qfree.chk subspace.qhom.chk			\
	subspace2 subspace2.qfree.chk subspace2.qhom.chk subspace2.sign

check-maxcutoff:
	@../check.template $(srcdir) "-maxcutoff" src/groebner/qsolve -q -m -o maxcutoff

check-maxinter:
	@../check.template $(srcdir) "-maxinter" src/groebner/qsolve -q -s -o maxinter 

check-local: check-maxcutoff check-maxinter
	@../flavorsum "-maxcutoff" "-maxinter"
	@../subdirsum $(SUBDIRS)

CLEANFILES = .num_failures .num_failures-maxcutoff .num_failures-maxinter

