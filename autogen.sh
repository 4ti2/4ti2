#!/bin/sh
aclocal -I m4
autoheader
autoconf
libtoolize --force || glibtoolize --force
mkdir -p swig
automake --add-missing
(test -e swig/autogen.sh && cd swig && sh ./autogen.sh)
