EXTRA_DIST =											\
	44.mat 44.hil.chk									\
	55.mat 55.hil.chk									\
	magic33.mat magic33.hil.chk								\
	reduce_gaussian-fedora-bug-911437.lat							\
	33.hil.chk 33.mat 33.rel 33.sign 33.zfree.chk						\
	a1.hil.chk a1.mat a1.rel a1.sign							\
	m33.hil.chk m33.mat m33.rel m33.sign							\
	non-pointed-2.hil.chk non-pointed-2.mat non-pointed-2.sign non-pointed-2.zfree.chk	\
	non-pointed.hil.chk non-pointed.mat non-pointed.sign non-pointed.zfree.chk		\
	zero-rows.hil.dont-chk zero-rows.mat zero-rows.sign zero-rows.zfree.dont-chk

check-local:
	@../check.template $(srcdir) "" src/zsolve/hilbert -q

CLEANFILES = .num_failures
