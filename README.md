4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 1998, 2002, 2006, 2015, 2025 4ti2 team.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

COMPILING 4ti2
==============

If you checked out the sources using git, first, to generate the build system,
you need Automake, Autoconf, and Libtool. Run the following command:

	./autogen.sh

Run the following commands with the 4ti2 directory:

	./configure --prefix=INSTALLATION-DIRECTORY
	make
	make check
    make install-exec

The final command will install 4ti2 in a directory tree below the
INSTALLATION-DIRECTORY that you gave with the first command.  If you
omit the --prefix option, `make install' will install 4ti2 in the
/usr/local hierarchy.

You will need glpk and gmp installed first (see below).

The first command, 'make', compiles all the executables. The second
command, 'make check', runs a lot of automatic checks. This will take
a while.  If a check fails, then please notify the 4ti2 team.

You will need gcc version 3.4 or higher.

You will need an installed version of glpk (linear programming software). See
the website http://www.gnu.org/software/glpk for more information. The
version 4.7 has been tested. If you do not have glpk installed or 4ti2 cannot
find glpk, then the compilation will fail saying that it cannot find the file
"glpk.h".  If you have installed glpk but not in a location that 4ti2 finds by
default, then you will need to invoke 

     ./configure --with-glpk=/ROOT/OF/GLPK/INSTALLATION/HIERARCHY

You will also need an installed version of gmp, The GNU MP
Bignum Library, with c++ support enabled (see http://www.swox.com/gmp/ for more
details).  Versions 4.2.1 and 4.1.4 have been tested. If you are compiling a
version of gmp from the source, make sure that you enable c++ support
(--enable-cxx configure option).  If you have
installed gmp but not in a location that 4ti2 finds by default, then you
will need to invoke

     ./configure --with-gmp=/ROOT/OF/GMP/INSTALLATION/HIERARCHY

If you have gmp but not with c++ support, then ./configure will
fail with an error saying that the file "gmpxx.h" cannot be found.


USING MACPORTS ON MAC OS X
==========================

Use the following commands.

     sudo port install gmp glpk
     ./configure --with-gmp=/opt/local --with-glpk=/opt/local
     make
     sudo make install 


INSTALLATION ON WINDOWS USING CYGWIN
====================================

1. Install Cygwin from https://www.cygwin.com/

   In the installer, select the following packages:

   Devel: gcc-core gcc-g++ make
   Math: glpk gmp libglpk-devel libgmp-devel

2. Make sure you unpack the 4ti2 sources into a
   C:\DIRECTORY\WITHOUT\SPACES\IN\IT
   (for example, C:\4ti2)

3. Open the Cygwin terminal

4. Type:

		cd /cygdrive/c/DIRECTORY/WITHOUT/SPACES/IN/IT
		./configure
		make
		make install

5. Now you can run 4ti2's commands from the Cygwin terminal.


DOCUMENTATION
=============

See the manual or the website https://4ti2.github.io/ for information on using 4ti2.


PYTHON INTERFACES
=================

- https://github.com/alfsan/Py4ti2

- https://pypi.org/project/passagemath-latte-4ti2/
