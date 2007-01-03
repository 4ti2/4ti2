#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(4ti2, 1.3)
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile src/groebner/Makefile.inc src/util/Makefile.inc 
			  src/ppi/Makefile.inc src/zsolve/Makefile
			  m4/Makefile])
AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX

## Check for specific compiler flags
AC_LANG(C++)
# Pick good default optimization flags
AX_CXX_MAXOPT
# Turn on most compiler warnings
AX_CXXFLAGS_WARN_ALL
# The flag -ftrapv means that arithmetic overflow checking is enabled.
# This only works with the g++ compiler version 3.4 and newer versions.
AX_CHECK_COMPILER_FLAGS(-ftrapv, TRAPV_FLAG=-ftrapv)
AC_SUBST(TRAPV_FLAG, ${TRAPV_FLAG})

AC_LANG(C)
AX_CC_MAXOPT
AX_CFLAGS_WARN_ALL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h sys/param.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

# Check for GNU MP.
AC_LANG(C++)
LB_CHECK_GMP(,,[AC_MSG_ERROR([GNU MP not found!])])
if test -z "GMP_HAVE_CXX" ; then 
  AC_MSG_ERROR([Unable to use the GNU MP C++ library.  
Please make sure that GNU MP was configured with "--enable-cxx".])
fi

LB_CHECK_GLPK(,,[AC_MSG_ERROR([GNU Linear Programming Kit not found!])])

AC_OUTPUT
