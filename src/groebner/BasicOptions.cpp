/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "BasicOptions.h"
#include "Globals.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

using namespace _4ti2_;

BasicOptions* BasicOptions::o = new BasicOptions;

BasicOptions::BasicOptions()
{
    output = VERBOSE;
}

BasicOptions*
BasicOptions::instance()
{
    return o;
}

void
BasicOptions::process_options(int argc, char** argv)
{
    int c;
    while (1) {
#ifdef _GNU_SOURCE
        int option_index = 0;
        static struct option long_options[] = {
            {"precision",        1, 0,'p'},
            {"quiet",            0, 0,'q'},
            {"help",             0, 0,'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "p:qh", long_options, &option_index);
#else
        c = getopt(argc, argv, "p:qh");
#endif
        if (c == -1)
            break;

        switch (c) {
        case 'p': // The precision (i.e. 32, 64, or arbitrary)
            if (std::string("32").find(optarg) == 0) { }
            else if (std::string("64").find(optarg) == 0) { }
            else if (std::string("arbitrary").find(optarg) == 0) { }
            else { unrecognised_option_argument("-p, --precision"); }
            break;
        case 'q':
            output = SILENT;
            out = new std::ofstream;
            break;
        case 'h':
        case '?':
        case ':':
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
BasicOptions::print_usage()
{
    std::cerr << "Usage: " << Globals::exec << " [options] <filename>\n\n";
    std::cerr << "\
Options:\n\
  -p, --precision=PREC       Select PREC as the integer arithmetic precision.\n\
                             PREC is one of the following: `64' (default),\n\
                             `32', and `arbitrary' (only `arb` is needed).\n\
  -q, --quiet                Do not output anything to the screen.\n\
  -h, --help                 Display this help and exit.\n\
\n";
}

void
BasicOptions::unrecognised_option_argument(const char* option)
{
   std::cerr << "4ti2: ";
   std::cerr << "Unrecognised argument \"" << optarg << "\" ";
   std::cerr << "for the option " << option << ".\n\n";
   print_usage();
   exit(1);
}
