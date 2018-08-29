EXTRA_DIST =					\
	33.mat 33.cir.chk			\
	33.1.mat 33.1.cir.chk 33.1.sign		\
	33.2.mat 33.2.cir.chk 33.2.sign		\
	33.3.mat 33.3.cir.chk 33.3.sign		\
	44.mat 44.cir.chk			\
	ppi3.mat ppi3.cir.chk

check-local:
	@../check.template $(srcdir) "" src/groebner/circuits -q -m
	@../check.template $(srcdir) "" src/groebner/circuits -q -s

CLEANFILES = .num_failures
