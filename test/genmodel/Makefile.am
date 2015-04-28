EXTRA_DIST =					\
	333.mod 333.mat.chk			\
	3x6.mod 3x6.mat.chk			\
	3x4x6.mod 3x4x6.mat.chk			\
	k2_3.mod k2_3.mat.chk

check-local:
	@../check.template $(srcdir) "" src/util/genmodel -q

CLEANFILES = .num_failures
