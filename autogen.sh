#!/bin/sh
aclocal -I m4
autoconf
libtoolize
automake --add-missing
