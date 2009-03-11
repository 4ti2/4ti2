AC_DEFUN([CHECK_TRAPV], 
[
# The flag -ftrapv means that arithmetic overflow checking is enabled.
# This only works with the g++ compiler version 3.4 and newer versions.
AX_CHECK_COMPILER_FLAGS(-ftrapv, TRAPV_FLAG=-ftrapv)
AC_SUBST(TRAPV_FLAG, ${TRAPV_FLAG})
# Check whether -ftrapv appears to actually trigger traps on overflows.
save_CXXFLAGS=${CXXFLAGS}
CXXFLAGS="${CXXFLAGS} ${TRAPV_FLAG}"
AC_MSG_CHECKING([whether -ftrapv actually seems to work for int])
trapv_int=no
AC_TRY_RUN([
#include<stdio.h>
volatile int a, b, c;	
int main()
{
  a = 1 << (sizeof(int) * 8 - 2);
  b = 1 << (sizeof(int) * 8 - 2);
  c = a + b;
  return 0;
}
], [trapv_int=no], [trapv_int=yes])
if test "$trapv_int" = "yes"; then
   AC_DEFINE(HAVE_TRAPV_INT, 1, [Define if -ftrapv is working for int])
fi
AC_MSG_RESULT($trapv_int)
AC_MSG_CHECKING([whether -ftrapv actually seems to work for long long])
trapv_long_long=no
AC_TRY_RUN([
volatile long long a, b, c;	
int main()
{
  a = 1 << (sizeof(long long) * 8 - 2);
  b = 1 << (sizeof(long long) * 8 - 2);
  c = a + b;
  return 0;
}
], [trapv_long_long=no], [trapv_long_long=yes])
if test "$trapv_long_long" = "yes"; then
   AC_DEFINE(HAVE_TRAPV_LONG_LONG, 1, [Define if -ftrapv is working for long long])
fi
AC_MSG_RESULT($trapv_long_long)
CXXFLAGS=${save_CXXFLAGS}

])
