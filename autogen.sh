#!/bin/sh
## gnulib-tool --import --lib=libgnu --source-base=lib --m4-base=m4 --doc-base=doc --tests-base=tests --aux-dir=. --no-conditional-dependencies --no-libtool --macro-prefix=gl vasprintf getopt-gnu
gnulib-tool --update
aclocal -I m4
autoheader
autoconf
libtoolize --force || glibtoolize --force
automake --add-missing
