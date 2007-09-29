EXTRA_DIST = \
	33.mat 33.zfree.chk 33.zhom.chk 33.zinhom.chk \
	hppi6.mat hppi6.lb hppi6.rel hppi6.sign hppi6.ub hppi6.zhom.chk hppi6.zinhom.chk \
    m33.mat m33.rel m33.sign 


check-local:
	@../check.template $(srcdir) "" src/zsolve/zsolve -q

CLEANFILES = .num_failures
