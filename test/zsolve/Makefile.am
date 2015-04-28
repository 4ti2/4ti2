SUBDIRS = api

EXTRA_DIST = \
	33.mat 33.zfree.chk 33.zhom.chk 33.zinhom.chk \
	hppi6.mat hppi6.lb hppi6.rel hppi6.sign hppi6.ub hppi6.zhom.chk hppi6.zinhom.chk \
	m33.mat m33.rel m33.sign \
	lattice_generation_bug.mat lattice_generation_bug.rel lattice_generation_bug.sign lattice_generation_bug.rhs lattice_generation_bug.zhom.chk lattice_generation_bug.zinhom.chk \
	a1.mat a1.rel a1.rhs a1.sign a1.zhom.chk a1.zinhom.chk \
	dutour-testcase-2013-08-21.mat dutour-testcase-2013-08-21.rel dutour-testcase-2013-08-21.rhs dutour-testcase-2013-08-21.sign dutour-testcase-2013-08-21.zhom.chk dutour-testcase-2013-08-21.zinhom.chk \
	non-pointed-2.mat non-pointed-2.sign non-pointed-2.zhom.chk																	       \
	non-pointed.mat non-pointed.sign non-pointed.zfree.chk non-pointed.zhom.chk															       \
	zero-rows.mat zero-rows.sign zero-rows.zfree.dont-chk zero-rows.zhom.dont-chk \
	affine.lat affine.sign affine.zhom.chk affine.zinhom.chk \
	system.mat system.rel system.sign system.zhom.chk system.zinhom.chk


check-local:
	@../check.template $(srcdir) "" src/zsolve/zsolve -q

CLEANFILES = .num_failures
