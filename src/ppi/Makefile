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

# Include makefile variable settings.
include Makefile.inc

# The source files to compile
SRCFILES = $(wildcard *.cpp)
# The object files
OBJFILES = $(SRCFILES:%.cpp=$(OBJDIR)/%.o)
# The dependency files.
DEPFILES = $(OBJFILES:%.o=%.d)
# The directory for the object and dependency files.
OBJDIR = ../../obj/ppi
# The executable to build.
EXECUTABLE = $(BINDIR)/ppi
# The binary directory where all executables are placed.
BINDIR = ../../bin

# The default rule.
default: $(EXECUTABLE) $(SCRIPTS)

# The target to compile the executable.
$(EXECUTABLE): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(COMPILER) $(CXXFLAGS) $(INCLUDE) $^ $(LIBINCLUDE) $(LIBS) -o $@
	@strip $(EXECUTABLE)

# The target to build object files from source files.
$(OBJDIR)/%.o: %.cpp
	$(COMPILER) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# The target to build dependency files from source files.
$(OBJDIR)/%.d: %.cpp
	@mkdir -p $(OBJDIR)
	$(COMPILER) -MM -MT'$@ $(@:.d=.o)' $(CXXFLAGS) $(INCLUDE) $< > $@

# The target to remove all automatically generated files.
clean:
	rm -rf $(OBJDIR) $(EXECUTABLE)

# Include all the dependency information from the dependency files.
-include $(DEPFILES)
