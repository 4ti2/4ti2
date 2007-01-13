#!/bin/sh
aclocal -I m4
autoconf
libtoolize --force
automake --add-missing
