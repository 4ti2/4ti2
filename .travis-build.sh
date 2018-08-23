#! /bin/bash
set -e
autoreconf -fi
case x$DOCKER in
    x)
	mkdir -p build
	cd build
	../configure --srcdir=.. && CAT_CHECKDIR_ON_ERROR=t make distcheck
	;;
    *i386*)
	docker run -i -v "${PWD}:/src" $DOCKER /bin/bash -c "linux32 --32bit i386 /src/.docker-build.sh"
	;;
    *)
	docker run -i -v "${PWD}:/src" $DOCKER /bin/bash -c "/src/.docker-build.sh"
	;;
esac
