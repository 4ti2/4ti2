EXTRA_DIST =								\
	33.mat 33.ray.chk						\
	44.mat 44.ray.chk						\
	55.mat 55.ray.chk						\
	33.r 33.r.ray.chk 33.r.sign					\
	44.r.mat 44.r.ray.chk 44.r.sign					\
	55.r 55.r.ray.chk 55.r.sign					\
	ex ex.ray.chk ex.rel ex.sign					\
	ppp.r ppp.r.qfree.chk ppp.r.ray.chk ppp.r.sign			\
	small small.ray.chk small.rel small.sign			\
	subspace subspace.qfree.chk subspace.ray.chk subspace.sign	\
	subspace2 subspace2.qfree.chk subspace2.ray.chk subspace2.sign

check-m:
	@../check.template $(srcdir) "-m" src/groebner/rays -q -m 

check-s:
	@../check.template $(srcdir) "-s" src/groebner/rays -q -s 

check-local: check-m check-s
	@../flavorsum "-m" "-s"

CLEANFILES = .num_failures .num_failures-m .num_failures-s
