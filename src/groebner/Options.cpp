/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

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
*/

#include "Options.h"
#include "Globals.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

using namespace _4ti2_;

Options* Options::o = new Options;

Options::Options()
{
    output = VERBOSE;
}

Options*
Options::instance()
{
    return o;
}

void
Options::process_options(int argc, char** argv)
{
    int c;
    optind = 1;
    while (1) {
#ifdef _GNU_SOURCE
        int option_index = 0;
        static struct option long_options[] = {
            {"algorithm",        1, 0,'a'},
            {"generation",       1, 0,'g'},
            {"minimal",          1, 0,'m'},
            {"auto-reduce-freq", 1, 0,'r'},
            {"output-freq",      1, 0,'f'},
            {"truncation",       1, 0,'t'},
            {"precision",        1, 0,'p'},
            {"quiet",            0, 0,'q'},
            {"help",             0, 0,'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "g:a:m:r:f:t:p:qh",
                 long_options, &option_index);
#else
        c = getopt(argc, argv, "g:a:m:r:f:t:p:qh");
#endif
        if (c == -1)
            break;

        switch (c) {
        case 'a':
            if (std::string("fifo").find(optarg) == 0)
            { Globals::algorithm = Globals::NORMAL; }
            else if (std::string("weighted").find(optarg) == 0)
            { Globals::algorithm = Globals::WEIGHTED; }
            else if (std::string("unbounded").find(optarg) == 0)
            { Globals::algorithm = Globals::GEBAUER_MOELLER; }
            else { unrecognised_option_argument("-a, --algorithm"); }
            break;
        case 'g':
            if (std::string("hybrid").find(optarg) == 0)
            { Globals::generation = Globals::HYBRID; }
            else if (std::string("project-and-lift").find(optarg) == 0)
            { Globals::generation = Globals::PROJECT_AND_LIFT; }
            else if (std::string("saturation").find(optarg) == 0)
            { Globals::generation = Globals::SATURATION; }
            else if (std::string("max-min").find(optarg) == 0)
            { Globals::generation = Globals::MAXMIN; }
            else { unrecognised_option_argument("-g, --generation"); }
            break;
        case 'm':
            if (std::string("yes").find(optarg) == 0)
            { Globals::minimal = true; }
            else if (std::string("no").find(optarg) == 0)
            { Globals::minimal = false; }
            else { unrecognised_option_argument("-m, --minimal"); }
            break;
        case 'r':
            if (sscanf(optarg, "%d", &Globals::auto_reduce_freq) != 1)
            {  unrecognised_option_argument("-r, --auto_reduce_freq"); }
            break;
        case 'f':
            if (sscanf(optarg, "%d", &Globals::output_freq) != 1)
            {  unrecognised_option_argument("-f, --output_freq"); }
            break;
        case 'q':
            output = SILENT;
            out = new std::ofstream;
            err = new std::ofstream;
            break;
        case 't':
            if (std::string("ip").find(optarg) == 0)
            { Globals::truncation = Globals::IP; }
            else if (std::string("lp").find(optarg) == 0)
            { Globals::truncation = Globals::LP; }
            else if (std::string("weight").find(optarg) == 0)
            { Globals::truncation = Globals::WEIGHT; }
            else if (std::string("none").find(optarg) == 0)
            { Globals::truncation = Globals::NONE; }
            else { unrecognised_option_argument("-t, --truncation"); }
            break;
        case 'p': // The precision (i.e. int32, int64, or arbitrary)
            if (std::string("32").find(optarg) == 0) { }
            else if (std::string("64").find(optarg) == 0) { }
            else if (std::string("arbitrary").find(optarg) == 0) { }
            else { unrecognised_option_argument("-p, --precision"); }
            break;
        case 'h': case '?': case ':':
            print_usage();
            exit(1);
            break;

        default:
            std::cerr << "Error: getopt returned unknown character code\n";
            print_usage();
            exit(1);
        }
    }

    if (optind == argc-1)
    {
        filename = argv[optind];
    }
    else
    {
        std::cerr << "Command Line Error: Incorrect number of arguments.\n";
        print_usage();
        exit(1);
    }
}

void
Options::print_usage()
{
    if (Globals::exec == "groebner")
    {
        std::cerr << "Usage: groebner [options] <PROJECT>\n\n";
        std::cerr << "Computes the Groebner basis of a lattice.\n";
        std::cerr << "\
Input Files:\n\
  PROJECT.mat         A matrix (optional if lattice basis is given).\n\
  PROJECT.lat         A lattice basis (optional if matrix is given).\n\
  PROJECT.cost        The cost matrix (optional, default is degrevlex).\n\
                      Ties are broken with degrevlex.\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative and '0' means a free variable).\n\
                      It is optional, and the default is all non-negative.\n\
  PROJECT.mar         The Markov basis/generating set of the lattice (optional).\n\
  PROJECT.weights     The weight vectors used for truncation (optional).\n\
  PROJECT.weights.max The maximum weights used for truncation.\n\
                      This file is needed when PROJECT.weights exists.\n\
  PROJECT.zsol        An integer solution to specify a fiber (optional).\n\
                      The integer solution is used for truncation.\n\
Output Files:\n\
  PROJECT.gro         The Groebner basis of the lattice.\n\n";
    }
    else if (Globals::exec == "markov")
    {
        std::cerr << "Computes the markov basis/generating set of a lattice.\n";
        std::cerr << "Usage: markov [options] <PROJECT>\n\n";
        std::cerr << "\
Input Files:\n\
  PROJECT             A matrix (optional only if lattice basis is given).\n\
  PROJECT.lat         A lattice basis (optional only if matrix is given).\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative and '0' means a free variable).\n\
                      It is optional, and the default is all non-negative.\n\
  PROJECT.weights     The weight vectors used for truncation (optional).\n\
  PROJECT.weights.max The maximum weights used for truncation.\n\
                      This file is needed when PROJECT.weights exists.\n\
  PROJECT.zsol        An integer solution to specify a fiber (optional).\n\
                      The integer solution is used for truncation.\n\
Output Files:\n\
  PROJECT.mar         The Markov basis/generating set of the lattice.\n";
    }
    else
    {
        std::cerr << "Usage: " << Globals::exec << " [options] <filename>\n\n";
    }
    std::cerr << "\
Options:\n\
  -p, --precision=PREC       Select PREC as the integer arithmetic precision.\n\
                             PREC is one of the following: `64' (default),\n\
                             `32', and `arbitrary' (only `arb` is needed).\n\
  -a, --algorithm=ALG        Select ALG as the completion procedure for\n\
                             computing Groebner bases. ALG is one of\n\
                             `fifo', `weighted', or 'unbounded.'\n\
  -g, --generation=ALG       Select ALG as the procedure for computing \n\
                             a generating set or Markov basis. ALG is\n\
                             one of `hybrid' (default), `project-and-lift',\n\
                             `max-min', or 'saturation'.\n\
  -t, --truncation=TRUNC     Set TRUNC as the truncation method.  TRUNC is\n\
                             of the following: `ip', `lp', `weight' (default),\n\
                             or `none'. Only relevant if `zsol' is given.\n\
  -m, --minimal=STATE        If STATE is `yes' (default), then 4ti2 will\n\
                             compute a minimal Markov basis. If STATE is\n\
                             'no', then the Markov basis will not \n\
                             necessarily be minimal.\n\
  -r, --auto-reduce-freq=n   Set the frequency of auto reduction.\n\
                             (default is 2500).\n\
  -f, --output-freq=n        Set the frequency of output (default is 1000).\n\
  -q, --quiet                Do not output anything to the screen.\n\
  -h, --help                 Display this help and exit.\n\
\n\
Only short options are supported on sun machines.\n\
\n";
}

void
Options::unrecognised_option_argument(const char* option)
{
   std::cerr << "4ti2: ";
   std::cerr << "Unrecognised argument \"" << optarg << "\" ";
   std::cerr << "for the option " << option << ".\n\n";
   print_usage();
   exit(1);
}
