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
#--# This Makefile has only been tested with GNU make.
#--#
#--# This Makefile requires perl perl, v5.6.1 to be installed. See
#--# the Makefile variable PERL which points to the executable.
#------------------------------------------------------------------
#--# Type
#--# 	make
#--# in order to see the available commands of this Makefile and a
#--# short description.
#------------------------------------------------------------------

# Convention for this file: top-level rules come before the rules that
# are used, so that one can basically read the Makefile from begin to
# end.

# Where are we? The make command should be started from the 
# top-level 4ti2 directory.
ROOT4ti2 := $(shell pwd)

include Makefile.inc

#------------------------------------------------------------------

# Some flags for the C compiler to produce optimized code.
CFLAGS = -O3 -Wall
##??CFLAGS2 = -g

# Additional flags for the linker.
LIBS = -lm -L. -l4ti2

COMPILER=$(CC) $(CFLAGS)

# Name of the library of all .o files.
LIB4ti2=lib4ti2.a

#------------------------------------------------------------------
#-- The commands ...
#------------------------------------------------------------------

#-- Compile everything
all: exe
#all: archive

exe: $(EXECUTABLES) install

#-- Generate executable files
$(EXECUTABLES): %: %.d %_main.c $(LIB4ti2)
	@echo Compiling executable $(*F) ...
	@$(COMPILER) $@_main.c -o $@ $(LIBS)

# Copies the executables into the bin directory.
install: $(EXECUTABLES)
	@mkdir -p $(BINDIR)
	cp -r $(EXECUTABLES) $(BINDIR)

#-- Create the library:
$(LIB4ti2): $(SRC:%=$(LIB4ti2)(%.o))


#-- Compile a file and put the .o file into the library.
$(SRC:%=$(LIB4ti2)(%.o)): $(LIB4ti2)(%.o): %.c
	@$(COMPILER) -c $<
	@$(AR) $(ARREPLFLAGS) $(LIB4ti2) $*.o
	@$(RM) $*.o


#-- Generate header dependencies and include them into the Makefile.
$(SRC:%=%.d): %.d: %.c
	@echo Generate header dependencies for $*.c ...
	$(CC) -MM -MT'$(LIB4ti2)($*.o) $@' $(CPPFLAGS) $< > $@


#-- This is a dummy generation of an empty file needed as a workaround
#-- for the following `include'.
Makefile.d:
	touch $@

include Makefile.d $(wildcard *.d)

# For each executable we generate a wrapper file that just contains a
# main function which calles the appropriate main function from the
# library, namely NAME_main().
# If there is no header file NAME.h we generate even that.
$(EXECUTABLES:%=%_main.c):
	@echo Generating source for $(@F) ...
	@echo "#include \"myheader.h\""			>  $@
	@echo "#include \"$(*F:%_main=%).h\""		>> $@
	@echo "int main(int argc, char *argv[]) {"	>> $@
	@echo "  return ($(*F)(argc, argv));"		>> $@
	@echo "}"					>> $@
	@(if [ ! -f $(*F).h ]; then \
		echo Generating $(*F).h ...; \
		echo "int $(*F)(int, char**);" > $(*F:%_main=%).h;\
	fi)

clean:
	rm -f $(SRC:%=%.d) $(SRC:%=%.o) $(LIB4ti2) $(EXECUTABLES) $(EXECUTABLES:%=%_main.c)

###################################################################
#EMACS Local Variables:
#EMACS mode: makefile
#EMACS End:
