# Check for GLPK, derived from gmp-check.m4, which comes from linbox

dnl LB_CHECK_GLPK ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for GLPK and define GLPK_CFLAGS and GLPK_LIBS

AC_DEFUN([LB_CHECK_GLPK],
[
DEFAULT_CHECKING_PATH="DEFAULT"

GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"

AC_ARG_WITH(glpk,
		AS_HELP_STRING([--with-glpk=DIR],
			[Use the GLPK library installed in DIR.
			 Otherwise, the library is searched in the standard locations 
			 (like "/usr" or "/usr/local").]),
		[if test "$withval" = yes ; then
			GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	         elif test "$withval" != no ; then
			GLPK_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	        fi],
		[GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}
BACKUP_CXX=${CXX}

# Following does not work because libtool is only created when configure has completed.
##CXX="./libtool --mode=link --tag=CXX ${CXX}"

AC_MSG_CHECKING(for GLPK)

for GLPK_HOME in ${GLPK_HOME_PATH} 
  do	
	if test "x$GLPK_HOME" == "xDEFAULT" -o -r "$GLPK_HOME/include/glpk.h"; then

		if test "x$GLPK_HOME" != "xDEFAULT" ; then
			GLPK_CFLAGS="-I${GLPK_HOME}/include"
			# Use this version during actual build:
			GLPK_LIBS="-L${GLPK_HOME}/lib -lglpk"
			### We used to use -R here: 
			### but we don't use this anywhere else;
			### it's probably obsolete
			### -R${GLPK_HOME}/lib 
			# During configure, we don't use libtool, 
			# so cannot portably use the -R option. 
			GLPK_LIBS_NOLIBTOOL="-L${GLPK_HOME}/lib -lglpk"
		else
			GLPK_CFLAGS=
			GLPK_LIBS="-lglpk"		
			GLPK_LIBS_NOLIBTOOL="-lglpk"		
		fi

		CXXFLAGS="${CXXFLAGS} ${GLPK_CFLAGS}"
		LIBS="${LIBS} ${GLPK_LIBS_NOLIBTOOL}"

		AC_LINK_IFELSE(AC_LANG_PROGRAM([extern "C" {
#include <glpk.h>
}], [LPX *lpx = lpx_create_prob(); lpx_delete_prob(lpx); ]),
		[
				AC_MSG_RESULT(found)
				AC_SUBST(GLPK_CFLAGS)
		  		AC_SUBST(GLPK_LIBS)
				AC_DEFINE(HAVE_GLPK,1,[Define if GLPK is installed])
				glpk_found="yes"
				ifelse([$2], , :, [$2])
				break
			],[			
				glpk_problem="$glpk_problem $GLPK_HOME"
				unset GLPK_CFLAGS
				unset GLPK_LIBS	
			])	
	fi
done

if test "x$glpk_found" != "xyes"; then
	AC_MSG_RESULT(not found)
	ifelse($3, , :, $3)
fi

CXX=${BACKUP_CXX}
CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
