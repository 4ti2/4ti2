#!/bin/sh
aclocal -I m4
autoheader
autoconf
libtoolize --force
automake --add-missing
(cd swig && sh ./autogen.sh)
