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

# This Makefile has only been tested with GNU make.

# Convention for this file: top-level rules come before the rules that
# are used, so that one can basically read the Makefile from beginning to
# end.

include Makefile.inc

#------------------------------------------------------------------

# Some flags for the C compiler to produce optimized code.
CFLAGS = -O3 -Wall
##??CFLAGS2 = -g

# Additional flags for the linker.
LIBS =

COMPILER=$(CC) $(CFLAGS)

#------------------------------------------------------------------
#-- The commands ...
#------------------------------------------------------------------

#-- Compile everything
all: exe
#all: archive

exe: $(EXECUTABLES)

#-- Generate executable files
$(EXECUTABLES): $(BINDIR)/%: $(OBJDIR)/%_main.o $(DEP) $(OBJ)
	@mkdir -p $(BINDIR)
	$(COMPILER) $< $(OBJ) $(LIBS) -o $@

#-- Compile a source file.
$(OBJDIR)/%.o: %.c
	@mkdir -p $(OBJDIR)
	$(COMPILER) -c $< -o $@

#-- Generate header dependencies and include them into the Makefile.
$(OBJDIR)/%.d: %.c
	@mkdir -p $(OBJDIR)
	@echo "Generating dependencies for $<"
	@$(CC) -MM -MT'$@ $(@:.d=.o)' $(CPPFLAGS) $< > $@

#-- This is a dummy generation of an empty file needed as a workaround
#-- for the following `include'.
Makefile.d:
	touch $@

include Makefile.d $(wildcard $(OBJDIR)/*.d)

# For each executable we generate a wrapper file that just contains a
# main function which calles the appropriate main function from the
# library, namely NAME_main().
# If there is no header file NAME.h we generate even that.
$(MAINSRC):
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
	rm -rf $(OBJDIR) $(EXECUTABLES) $(MAINSRC)

###################################################################
#EMACS Local Variables:
#EMACS mode: makefile
#EMACS End:
