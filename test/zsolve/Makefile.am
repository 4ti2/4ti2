EXTRA_DIST = \
	33.mat 33.zfree.chk 33.zhom.chk 33.zinhom.chk \
	hppi6.mat hppi6.lb hppi6.rel hppi6.sign hppi6.ub hppi6.zhom.chk hppi6.zinhom.chk \
	m33.mat m33.rel m33.sign \
	lattice_generation_bug.mat lattice_generation_bug.rel lattice_generation_bug.sign lattice_generation_bug.lat lattice_generation_bug.rhs lattice_generation_bug.zhom.chk lattice_generation_bug.zinhom.chk \
	a1.mat a1.rel a1.rhs a1.sign a1.zhom.chk a1.zinhom.chk

check-local:
	@../check.template $(srcdir) "" src/zsolve/zsolve -q

CLEANFILES = .num_failures
