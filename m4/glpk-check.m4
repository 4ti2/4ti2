# Check for GLPK, derived from gmp-check.m4, which comes from linbox

dnl LB_CHECK_GLPK ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GLPK_CFLAGS and GLPK_LIBS

AC_DEFUN([LB_CHECK_GLPK],
[
DEFAULT_CHECKING_PATH="/usr /usr/local"

GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"

AC_ARG_WITH(glpk,
		[  --with-glpk= <path>|yes|no
	   				   Use GLPK library. 
					   If argument is no, you do not have the library installed on your machine.
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path "/usr" or "/usr/local"  (set as default).
	 				   Otherwise you give the <path> to the directory which contain the library. 
		],
		[if test "$withval" = yes ; then
			GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	         elif test "$withval" != no ; then
			GLPK_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	        fi],
		[GLPK_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for GLPK)

for GLPK_HOME in ${GLPK_HOME_PATH} 
  do	
	if test -r "$GLPK_HOME/include/glpk.h"; then

		if test "x$GLPK_HOME" != "x/usr" -a "x$GLPK_HOME" != "x/usr/local"; then
			GLPK_CFLAGS="-I${GLPK_HOME}/include"
			GLPK_LIBS="-L${GLPK_HOME}/lib -R${GLPK_HOME}/lib -lglpk"	
		else
			GLPK_CFLAGS=
			GLPK_LIBS="-lglpk"		
		fi

		CXXFLAGS="${CXXFLAGS} ${GLPK_CFLAGS}"
		LIBS="${LIBS} ${GLPK_LIBS}"

		AC_TRY_LINK_FUNC(glp_lpx_create_prob,
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

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
