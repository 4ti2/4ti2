# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Peter Malkin.
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

bin_PROGRAMS = 4ti2int32 4ti2int64 
if HAVE_GMP_WITH_CXX
bin_PROGRAMS += 4ti2gmp
endif

lib_LTLIBRARIES = lib4ti2int32.la lib4ti2int64.la lib4ti2gmp.la

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(srcdir)/..

# 4ti2 uses the glpk linear programming library (GLPL).
# Add include path to the file glpk.h if necessary.
AM_CXXFLAGS += $(GLPK_CFLAGS)
# Add library include path for the library libglpk if necessary.
LIBS += $(GLPK_LIBS)

# 4ti2 needs the next option to access maximum values of integer types.
AM_CXXFLAGS += -D__STDC_LIMIT_MACROS
# The optimisation flags.
AM_CXXFLAGS += -DNDEBUG
# Uncomment the next line for debugging.
#AM_CXXFLAGS += -g
# Uncomment the next line for profiling.
#AM_CXXFLAGS += -pg

lib4ti2sources =				\
	Algorithm.cpp				\
	BasicCompletion.cpp			\
	BasicGeneration.cpp			\
	BasicOptions.cpp			\
	BasicReduction.cpp			\
	BinomialArray.cpp			\
	BinomialArrayStream.cpp			\
	BinomialCollection.cpp			\
	Binomial.cpp				\
	BinomialFactory.cpp			\
	BinomialSet.cpp				\
	BinomialSetStream.cpp			\
	BinomialStream.cpp			\
	BitSetStream.cpp			\
	Bounded.cpp				\
	CircuitOptions.cpp			\
	Completion.cpp				\
	DiagonalAlgorithm.cpp			\
	EliminateAlgorithm.cpp			\
	Euclidean.cpp				\
	Extract.cpp				\
	Feasible.cpp				\
	FeasibleStream.cpp			\
	FilterReduction.cpp			\
	FlipCompletion.cpp			\
	GeneratingSet.cpp			\
	Globals.cpp				\
	GroebnerBasis.cpp			\
	HermiteAlgorithm.cpp			\
	HybridGenSet.cpp			\
	LatticeBasis.cpp			\
	LongDenseIndexSet.cpp			\
	LongDenseIndexSetStream.cpp		\
	Markov.cpp				\
	MaxMinGenSet.cpp			\
	Minimize.cpp				\
	MinimizeOptions.cpp			\
	OnesReduction.cpp			\
	Optimise.cpp				\
	Options.cpp				\
	OrderedCompletion.cpp			\
	ProjectLiftGenSet.cpp			\
	QSolveAlgorithm.cpp 		\
	SaturationGenSet.cpp			\
	ShortDenseIndexSet.cpp			\
	ShortDenseIndexSetStream.cpp		\
	Statistics.cpp				\
	SyzygyCompletion.cpp			\
	SyzygyGeneration.cpp			\
	Timer.cpp				\
	Truncated.cpp				\
	VectorArray.cpp				\
	VectorArrayStream.cpp			\
	Vector.cpp				\
	VectorStream.cpp			\
	WalkAlgorithm.cpp			\
	WalkOptions.cpp				\
	WeightAlgorithm.cpp			\
	WeightedBinomialSet.cpp			\
	WeightedReduction.cpp

groebnerincludedir=$(includedir)/groebner

groebnerinclude_HEADERS =			\
	Algorithm.h				\
	BasicCompletion.h			\
	BasicGeneration.h			\
	BasicOptions.h				\
	BasicReduction.h			\
	BinomialArray.h				\
	BinomialArrayStream.h			\
	BinomialCollection.h			\
	BinomialFactory.h			\
	Binomial.h				\
	BinomialSet.h				\
	BinomialSetStream.h			\
	BinomialStream.h			\
	BitSet.h				\
	BitSetStream.h				\
	Bounded.h				\
	CircuitImplementation.h			\
	CircuitMatrixAlgorithm.h		\
	CircuitOptions.h			\
	CircuitSupportAlgorithm.h		\
	Completion.h				\
	DataType.h				\
	Debug.h					\
	DiagonalAlgorithm.h			\
	EliminateAlgorithm.h			\
	Euclidean.h				\
	Extract.h				\
	Feasible.h				\
	FeasibleStream.h			\
	Filter.h				\
	FilterReduction.h			\
	FlipCompletion.h			\
	GeneratingSet.h				\
	Generation.h				\
	Globals.h				\
	Grading.h				\
	GroebnerBasis.h				\
	HermiteAlgorithm.h			\
	HybridGenSet.h				\
	IndexBinomialSet.h			\
	Index.h					\
	IndexSetConverter.h			\
	LatticeBasis.h				\
	LongDenseIndexSet.h			\
	LongDenseIndexSetStream.h		\
	Markov.h				\
	MaxMinGenSet.h				\
	Minimize.h				\
	MinimizeOptions.h			\
	OnesReduction.h				\
	OnesTree.h				\
	Optimise.h				\
	Options.h				\
	OrderedCompletion.h			\
	Permutation.h				\
	ProjectLiftGenSet.h			\
	QSolveAlgorithm.h		\
	QSolveConsOrder.h		\
	QSolveVariant.h			\
	RayImplementation.h			\
	RayMatrixAlgorithm.h			\
	RaySupportAlgorithm.h			\
	Reduction.h				\
	SaturationGenSet.h			\
	ShortDenseIndexSet.h			\
	ShortDenseIndexSetStream.h		\
	SignPattern.h				\
	Size.h					\
	Statistics.h				\
	SupportTree.h				\
	SyzygyCompletion.h			\
	SyzygyGeneration.h			\
	TermOrder.h				\
	Timer.h					\
	Truncate.h				\
	VectorArray.h				\
	VectorArrayStream.h			\
	Vector.h				\
	VectorStream.h				\
	WalkAlgorithm.h				\
	WalkOptions.h				\
	WeightAlgorithm.h			\
	WeightedBinomialSet.h			\
	WeightedReduction.h			\
	Weight.h

# These files are C++ files #include'd in one of the above headers
groebnerinclude_HEADERS +=			\
	CircuitImplementation.tpp		\
	CircuitMatrixAlgorithm.tpp		\
	CircuitSupportAlgorithm.tpp		\
	DiagonalAlgorithm.tpp			\
	HermiteAlgorithm.tpp			\
	OnesTree.tpp				\
	RayImplementation.tpp			\
	RayMatrixAlgorithm.tpp			\
	RaySupportAlgorithm.tpp			\
	SupportTree.tpp

4ti2sources =					\
	circuits_main.cpp			\
	groebner_main.cpp			\
	markov_main.cpp				\
	minimize_main.cpp			\
	normalform_main.cpp			\
	qsolve_main.cpp				\
	rays_main.cpp				\
	walk_main.cpp				\
	zbasis_main.cpp				\
	main.cpp

# These headers are not installed
noinst_HEADERS = 				\
	circuits_main.h				\
	groebner_main.h				\
	markov_main.h				\
	minimize_main.h				\
	normalform_main.h			\
	qsolve_main.h				\
	rays_main.h				\
	walk_main.h				\
	zbasis_main.h

WRAPPERSCRIPTS =				\
	circuits				\
	groebner				\
	markov					\
	minimize				\
	normalform				\
	qsolve					\
	rays					\
	walk					\
	zbasis

nodist_bin_SCRIPTS = $(WRAPPERSCRIPTS)
DISTCLEANFILES = $(WRAPPERSCRIPTS)

# 16 bit precision flags.
# 4ti2int16_LDADD = lib4ti2int16.la
# 4ti2int16_CPPFLAGS = -D_4ti2_INT16_
# 4ti2int16_SOURCES = $(4ti2sources)
# lib4ti2int16_la_CPPFLAGS = -D_4ti2_INT16_
# lib4ti2int16_la_SOURCES = $(lib4ti2sources)

# 32 bit precision flags.
4ti2int32_LDADD = lib4ti2int32.la
4ti2int32_CPPFLAGS = -D_4ti2_INT32_
4ti2int32_CXXFLAGS = $(TRAPV) $(AM_CXXFLAGS)
4ti2int32_SOURCES = $(4ti2sources)
lib4ti2int32_la_CPPFLAGS = -D_4ti2_INT32_
lib4ti2int32_la_CXXFLAGS = $(TRAPV) $(AM_CXXFLAGS)
lib4ti2int32_la_SOURCES = $(lib4ti2sources)

# 64 bit precision flags.
4ti2int64_LDADD = lib4ti2int64.la
4ti2int64_CPPFLAGS = -D_4ti2_INT64_
4ti2int64_CXXFLAGS = $(TRAPV) $(AM_CXXFLAGS)
4ti2int64_SOURCES = $(4ti2sources)
lib4ti2int64_la_CPPFLAGS = -D_4ti2_INT64_
lib4ti2int64_la_CXXFLAGS = $(TRAPV) $(AM_CXXFLAGS)
lib4ti2int64_la_SOURCES = $(lib4ti2sources)

# Arbitrary precision flags.
# 4ti2 uses GMP (GLPL), an arbitrary precision arithmetic library.
if HAVE_GMP_WITH_CXX
4ti2gmp_LDADD = lib4ti2gmp.la
4ti2gmp_CPPFLAGS = -D_4ti2_GMP_ $(GMP_CFLAGS)
4ti2gmp_SOURCES = $(4ti2sources)
lib4ti2gmp_la_CPPFLAGS = -D_4ti2_GMP_
lib4ti2gmp_la_CPPFLAGS += $(GMP_CFLAGS)
lib4ti2gmp_la_LIBADD = $(GMP_LIBS)
lib4ti2gmp_la_SOURCES = $(lib4ti2sources)
endif

EXTRA_DIST = script.template

$(WRAPPERSCRIPTS): script.template
	cp $< $@

# Test and profiling
if HAVE_GMP_WITH_CXX
noinst_PROGRAMS = profile_qsolve
profile_qsolve_LDADD = lib4ti2gmp.la
profile_qsolve_CPPFLAGS = -D_4ti2_GMP_ $(GMP_CFLAGS)
profile_qsolve_SOURCES = profile_qsolve.cpp
endif
