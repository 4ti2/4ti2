# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Raymond Hemmecke, Ralf Hemmecke, Peter Malkin.
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
# The C compiler we use.
#CC=cc
CC=g++

#------------------------------------------------------------------

# The names of the final executables.
EXECNAMES= \
  genmodel \
  gensymm \
  output \

# The place to put the executables.
BINDIR = ../../bin
# The executables to create.
EXECUTABLES=$(EXECNAMES:%=$(BINDIR)/%)

# The source files.
SRC= \
  genmodel.c \
  gensymm.c \
  orbit.c \
  output.c \
  print.c \
  vector.c

# The main source files for each of the executables.
MAINSRC = $(EXECNAMES:%=%_main.c)

# The directory for the object and dependency files.
OBJDIR = ../../obj/util
# The object files to create.
OBJ = $(SRC:%.c=$(OBJDIR)/%.o)

# The dependency files to create.
DEP = $(OBJ:%.o=%.d) $(MAINSRC:%.c=$(OBJDIR)/%.d)

###################################################################
#EMACS Local Variables:
#EMACS mode: makefile
#EMACS End:
