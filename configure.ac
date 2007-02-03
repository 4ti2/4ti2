#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(4ti2, 1.3.1+autoconf)
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile src/Makefile 
		 src/groebner/Makefile src/util/Makefile 
		 src/ppi/Makefile src/zsolve/Makefile
		 m4/Makefile])
AC_CONFIG_FILES([test/checkall],		[chmod +x test/checkall])
AC_CONFIG_FILES([test/check.template], 		[chmod +x test/check.template])
AC_CONFIG_FILES([test/ppi/check-ppi], 		[chmod +x test/ppi/check-ppi])
AC_CONFIG_FILES([test/circuits/check], 		[chmod +x test/circuits/check])
AC_CONFIG_FILES([test/genmodel/check], 		[chmod +x test/genmodel/check])
AC_CONFIG_FILES([test/gensymm/check], 		[chmod +x test/gensymm/check])
AC_CONFIG_FILES([test/graver/check], 		[chmod +x test/graver/check])
AC_CONFIG_FILES([test/groebner/check], 		[chmod +x test/groebner/check])
AC_CONFIG_FILES([test/groebner/gmp/check], 	[chmod +x test/groebner/gmp/check])
AC_CONFIG_FILES([test/hilbert/check], 		[chmod +x test/hilbert/check])
AC_CONFIG_FILES([test/markov/check], 		[chmod +x test/markov/check])
AC_CONFIG_FILES([test/markov/gmp/check], 	[chmod +x test/markov/gmp/check])
AC_CONFIG_FILES([test/markov/saturation/check], [chmod +x test/markov/saturation/check])
AC_CONFIG_FILES([test/minimize/check], 		[chmod +x test/minimize/check])
AC_CONFIG_FILES([test/normalform/check],	[chmod +x test/normalform/check])
AC_CONFIG_FILES([test/output/binomials/check],	[chmod +x test/output/binomials/check])
AC_CONFIG_FILES([test/output/maple/check],	[chmod +x test/output/maple/check])
AC_CONFIG_FILES([test/output/positive/check],	[chmod +x test/output/positive/check])
AC_CONFIG_FILES([test/ppi/check],		[chmod +x test/ppi/check])
AC_CONFIG_FILES([test/qsolve/check],		[chmod +x test/qsolve/check])
AC_CONFIG_FILES([test/qsolve/gmp/check],	[chmod +x test/qsolve/gmp/check])
AC_CONFIG_FILES([test/rays/check], 		[chmod +x test/rays/check])
AC_CONFIG_FILES([test/walk/check], 		[chmod +x test/walk/check])
AC_CONFIG_FILES([test/zbasis/check], 		[chmod +x test/zbasis/check])
AC_CONFIG_FILES([test/zsolve/check], 		[chmod +x test/zsolve/check])

## Using the C compiler
AC_PROG_CC
AX_COMPILER_VENDOR
AX_CC_MAXOPT
AX_CFLAGS_WARN_ALL

## Using the C++ compiler
AC_PROG_CXX
# Check for specific compiler flags
AC_LANG(C++)
# Pick good default optimization flags
AX_COMPILER_VENDOR
AX_CXX_MAXOPT
# Turn on most compiler warnings
AX_CXXFLAGS_WARN_ALL
# The flag -ftrapv means that arithmetic overflow checking is enabled.
# This only works with the g++ compiler version 3.4 and newer versions.
AX_CHECK_COMPILER_FLAGS(-ftrapv, TRAPV_FLAG=-ftrapv)
AC_SUBST(TRAPV_FLAG, ${TRAPV_FLAG})

# Change the default behaviour to disable shared libraries.  The user
# may still override this default by specifying `--enable-shared'.
# The reasons for the changed default is that the complexities of
# building shared libraries are avoided, and also that static build
# can result in faster code on some architectures.
AC_DISABLE_SHARED
# This triggers selection of some optimization flags (-g -O2), 
# so it should appear _after_ the AX_CXX_MAXOPT, AX_CC_MAXOPT calls.
AC_PROG_LIBTOOL

# Checks for header files.
#AC_HEADER_STDC
#AC_CHECK_HEADERS([malloc.h stdlib.h string.h sys/param.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_CONST

# Checks for library functions.
#AC_FUNC_MALLOC

# Check for GNU MP.
AC_LANG(C++)
LB_CHECK_GMP(,,[AC_MSG_WARN([GNU MP not found.
4ti2 will be built without support for arbitrary-precision computations.])])
if test -z "${GMP_HAVE_CXX}" ; then 
  AC_MSG_WARN([Unable to use the GNU MP C++ library.  
Please make sure that GNU MP was configured with "--enable-cxx".
4ti2 will be built without support for arbitrary-precision computations.])
fi
AM_CONDITIONAL([HAVE_GMP_WITH_CXX], [test -n "${GMP_HAVE_CXX}"])

LB_CHECK_GLPK(,,[AC_MSG_ERROR([GNU Linear Programming Kit not found!])])

AC_OUTPUT
