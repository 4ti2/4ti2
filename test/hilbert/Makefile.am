EXTRA_DIST =					\
	44.mat 44.hil.chk			\
	55.mat 55.hil.chk			\
	magic33.mat magic33.hil.chk		\
	reduce_gaussian-fedora-bug-911437.lat

check-local:
	@../check.template $(srcdir) "" src/zsolve/hilbert -q

CLEANFILES = .num_failures
