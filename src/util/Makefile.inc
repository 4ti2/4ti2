# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 

#------------------------------------------------------------------
#--# This Makefile lists the names of files that should be
#--# compiled together with some common variables. It is meant to be
#--# included in the top Makefile.
#------------------------------------------------------------------

#------------------------------------------------------------------
#-- The programs that we use within the Makefile
#------------------------------------------------------------------
# Path to the perl executable.
PERL=perl

# The C compiler we use.
#CC=cc
CC=g++

# Create the directory including missing parent dirs.
MKDIR=mkdir

# Copy files without changing their timestamp.
CP=cp -p

# To create the release archive we need to put a script together with
# the tgz archive.
CAT=cat

# To produce the archive... (flags for replacing an archive member)
AR=ar
ARREPLFLAGS=rv

# Remove a file.
RM=rm

# Pack into an archive
TAR=tar
TARFLAGS=cvf

# Unpack
UNTAR=tar
UNTARFLAGS=xvf

# Compress the archive
GZIP=gzip
GZIPFLAGS=

# Uncompress the archive
GUNZIP=gunzip
GUNZIPFLAGS=

# Search regular expressions in a text file
GREP=grep

# Someone should know how to compile the documentation
LATEX=latex

#------------------------------------------------------------------

# Where are we? The make command should be started from the top-level
# 4ti2 directory.
ROOT4ti2 := $(shell pwd)

# The names of the final executables.
EXECUTABLES= \
  genmodel \
  gensymm \
  output \

# The source files that should go to the library.
SRC= \
  genmodel \
  gensymm \
  orbit \
  output \
  print \
  vector

# The place to put the executables.
BINDIR = ../../bin

###################################################################
#EMACS Local Variables:
#EMACS mode: makefile
#EMACS End:
