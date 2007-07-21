EXTRA_DIST =					\
	333 333.gra.chk				\
	bayer.gra.chk bayer.lat			\
	hppi4 hppi4.gra.chk			\
	hppi5 hppi5.gra.chk			\
	hppi6 hppi6.gra.chk			\
	hppi6.lb hppi6.rel hppi6.sign hppi6.ub	\
	small.gra.chk small.mat small.sign


check-local:
	@../check.template $(srcdir) "" src/zsolve/graver -q 

CLEANFILES = .num_failures
