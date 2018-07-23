#!/bin/sh
aclocal -I m4
autoheader
autoconf
libtoolize --force || glibtoolize --force
automake --add-missing
(cd swig && sh ./autogen.sh)
