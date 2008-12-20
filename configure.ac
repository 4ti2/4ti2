#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(4ti2, 1.3.1+autoconf-0.1)
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile src/Makefile src/4ti2/Makefile
		 src/groebner/Makefile src/util/Makefile 
		 src/ppi/Makefile src/zsolve/Makefile
		 m4/Makefile])

AC_CONFIG_HEADERS(src/config.h)

AC_CONFIG_FILES([test/Makefile])
AC_CONFIG_FILES([test/circuits/Makefile])
AC_CONFIG_FILES([test/genmodel/Makefile])
AC_CONFIG_FILES([test/gensymm/Makefile])
AC_CONFIG_FILES([test/graver/Makefile])
AC_CONFIG_FILES([test/groebner/Makefile])
AC_CONFIG_FILES([test/groebner/gmp/Makefile])
AC_CONFIG_FILES([test/hilbert/Makefile])
AC_CONFIG_FILES([test/markov/Makefile])
AC_CONFIG_FILES([test/markov/gmp/Makefile])
AC_CONFIG_FILES([test/markov/saturation/Makefile])
AC_CONFIG_FILES([test/markov/maxmin/Makefile])
AC_CONFIG_FILES([test/minimize/Makefile])
AC_CONFIG_FILES([test/normalform/Makefile])
AC_CONFIG_FILES([test/output/Makefile])
AC_CONFIG_FILES([test/output/binomials/Makefile])
AC_CONFIG_FILES([test/output/maple/Makefile])
AC_CONFIG_FILES([test/output/positive/Makefile])
AC_CONFIG_FILES([test/ppi/Makefile])
AC_CONFIG_FILES([test/qsolve/Makefile])
AC_CONFIG_FILES([test/qsolve/api/Makefile])
AC_CONFIG_FILES([test/qsolve/gmp/Makefile])
AC_CONFIG_FILES([test/rays/Makefile])
AC_CONFIG_FILES([test/walk/Makefile])
AC_CONFIG_FILES([test/zbasis/Makefile])
AC_CONFIG_FILES([test/zsolve/Makefile])
AC_CONFIG_FILES([test/zsolve/api/Makefile])

AC_CONFIG_FILES([test/check.template], 		[chmod +x test/check.template])
AC_CONFIG_FILES([test/subdirsum], 		[chmod +x test/subdirsum])
AC_CONFIG_FILES([test/flavorsum], 		[chmod +x test/flavorsum])
AC_CONFIG_FILES([test/reportfailures], 		[chmod +x test/reportfailures])
AC_CONFIG_FILES([test/ppi/check-ppi], 		[chmod +x test/ppi/check-ppi])

## Check what components to enable
AC_ARG_ENABLE(groebner, AS_HELP_STRING([--disable-groebner], [Disable building the groebner components]),
			[enable_groebner=${enableval}], [enable_groebner=yes])
AC_ARG_ENABLE(zsolve, AS_HELP_STRING([--disable-zsolve], [Disable building the zsolve components]),
			[enable_zsolve=${enableval}], [enable_zsolve=yes])
AM_CONDITIONAL([ENABLE_GROEBNER], [ test x${enable_groebner} == xyes ])
AM_CONDITIONAL([ENABLE_ZSOLVE], [ test x${enable_zsolve} == xyes ])

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
CHECK_TRAPV

### Taken from GMP's configure.in:
# The dead hand of AC_REQUIRE makes AC_PROG_LIBTOOL expand and execute
# AC_PROG_F77, even when F77 is not in the selected with_tags.  This is
# probably harmless, but it's unsightly and bloats our configure, so pretend
# AC_PROG_F77 has been expanded already.
#
# FIXME: Rumour has it libtool will one day provide a way for a configure.in
# to say what it wants from among supported languages etc.
#
AC_PROVIDE([AC_PROG_F77])

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

if test x${enable_groebner} == xyes ; then
# Check for GNU MP.
AC_LANG(C++)
LB_CHECK_GMP(,,[AC_MSG_WARN([GNU MP not found.
4ti2 will be built without support for arbitrary-precision computations.])])
if test -z "${GMP_HAVE_CXX}" ; then 
  AC_MSG_WARN([Unable to use the GNU MP C++ library.  
Please make sure that GNU MP was configured with "--enable-cxx".
4ti2 will be built without support for arbitrary-precision computations.])
fi

# Check for GLPK.
LB_CHECK_GLPK(,,[
  AC_MSG_ERROR([GNU Linear Programming Kit not found!
It is required for building the Groebner component of 4ti2.
To just build the zsolve component use --disable-groebner.])])

fi

AM_CONDITIONAL([HAVE_GMP_WITH_CXX], [test -n "${GMP_HAVE_CXX}"])

if test -n "${GMP_HAVE_CXX}" ; then 
  # Check whether we can convert between int32_t and int64_t and mpz_class (for the API).
  save_LIBS="${LIBS}"
  LIBS="${LIBS} ${GMP_LIBS}"
  AC_MSG_CHECKING(whether we can convert between int32_t and mpz_class)
  AC_TRY_LINK(
    [  #include <gmpxx.h>
    ],
    [  int32_t x = 815;
       mpz_class y(x);
    ],
    [ HAVE_MPZ_INT32_CONVERSION=yes;
      AC_DEFINE(HAVE_MPZ_INT32_CONVERSION, 1, [Define if we can convert between int32_t and mpz_class]) ],
    [ HAVE_MPZ_INT32_CONVERSION=no ])
  AC_MSG_RESULT(${HAVE_MPZ_INT32_CONVERSION})
  AC_MSG_CHECKING(whether we can convert between int64_t and mpz_class)
  AC_TRY_LINK(
    [  #include <gmpxx.h>
    ],
    [  int64_t x = 815;
       mpz_class y(x);
    ],
    [ HAVE_MPZ_INT64_CONVERSION=yes;
      AC_DEFINE(HAVE_MPZ_INT64_CONVERSION, 1, [Define if we can convert between int64_t and mpz_class]) ],
    [ HAVE_MPZ_INT64_CONVERSION=no ])
  AC_MSG_RESULT(${HAVE_MPZ_INT64_CONVERSION})
  LIBS="${save_LIBS}"
fi

## Check whether to build high-level-language interfaces using SWIG
AC_ARG_ENABLE(swig, AS_HELP_STRING([--enable-swig], [Enable building high-level-language interfaces using SWIG]),
		    [enable_swig=${enableval}], [enable_swig=no])
if test "$enable_swig" = yes; then
  # Call configure separately in this subdirectory.
  AC_CONFIG_SUBDIRS(swig)
  # If we use SWIG, we definitely want shared libraries.
  AC_ENABLE_SHARED
fi
AM_CONDITIONAL([ENABLE_SWIG], [ test x${enable_swig} == xyes ])

AC_OUTPUT
