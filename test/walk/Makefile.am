EXTRA_DIST =												\
	344.mat 344.cost 344.cost.start 344.gro.chk 344.gro.start 344.mar				\
	cuww1.mat cuww1.cost cuww1.cost.start cuww1.gro.chk cuww1.gro.start				\
	grin.mat grin.cost grin.cost.start grin.gro.chk grin.gro.start grin.mar				\
	prob20.mat prob20.cost prob20.cost.start prob20.gro.chk prob20.gro.start prob20.lat prob20.sign

check-32:
	@../check.template $(srcdir) "-32" src/groebner/walk -q 

check-arb:
	@../check.template $(srcdir) "-arb" src/groebner/walk -p arb -q 

check-local: check-32 check-arb
	@../flavorsum "-32" "-arb"

CLEANFILES = .num_failures .num_failures-32 .num_failures-arb
