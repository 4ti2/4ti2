AC_PREREQ(2.59)
AC_INIT(4ti2-swig, 1.3.1+autoconf-0.1)
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile allegrocommonlisp/Makefile])

# Only build shared libraries.
AC_DISABLE_STATIC
AC_ENABLE_SHARED

AC_PROG_CC
AC_PROG_CXX
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
AC_PROG_LIBTOOL

AC_ARG_WITH(swig,
    [  --with-swig[=ARG]       Specify location of SWIG],
    [SWIG=$withval],
    [AC_MSG_ERROR([Finding the right SWIG version is too hairy, please use --with-swig=ARG])])
AC_SUBST(SWIG)

AC_ARG_ENABLE(swig-allegrocommonlisp,
    [  --enable-swig-allegrocommonlisp    Enable build of Allegro Common Lisp FFI])
AM_CONDITIONAL(ENABLE_SWIG_ALLEGROCOMMONLISP, test x$enable_swig_allegrocommonlisp = xyes)

AC_OUTPUT
