EXTRA_DIST =						\
	333.mat 333.mar.chk				\
	334.mat 334.mar.chk				\
	cuww1.mat cuww1.bs cuww1.mar.chk cuww1.vs	\
	cuww2.mat cuww2.mar.chk				\
	cuww3.mat cuww3.A cuww3.mar.chk			\
	cuww4.mat cuww4.mar.chk				\
	cuww5.mat cuww5.mar.chk				\
	prob02.mat prob02.gro prob02.mar.chk		\
	prob04.mat prob04.gro prob04.mar.chk		\
	prob06.mat prob06.gro prob06.lat prob06.mar.chk	\
	prob08.mat prob08.gro prob08.mar.chk		\
	prob10.mat prob10.gro prob10.mar.chk		\
	prob12.mat prob12.gro prob12.mar.chk		\
	prob14.mat prob14.gro prob14.mar.chk		\
	prob16.mat prob16.gro prob16.lat prob16.mar.chk	\
	prob18.mat prob18.gro prob18.mar.chk		\
	prob20.mat prob20.mar.chk

check-local:
	@../../check.template $(srcdir) "" src/groebner/markov -q -gsaturation

CLEANFILES = .num_failures
