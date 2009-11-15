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

#include <iostream>
#include <fstream>
#include <unistd.h>
#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

#include "Options.h"
#include "Globals.h"

using namespace _4ti2_;

Options::Options(int argc, char** argv)
{
    optind = 1;
    precision = _4ti2_PREC_INT_64;
    int c;
    while (1) {
#ifdef _GNU_SOURCE
        int option_index = 0;
        static struct option long_options[] = {
            {"matrix",       0, 0,'m'},
            {"support",      0, 0,'s'},
            {"order",        1, 0,'o'},
            {"output-freq",  1, 0,'f'},
            {"precision",    1, 0,'p'},
            {"threads",      1, 0,'j'},
            {"quiet",        0, 0,'q'},
            {"help",         0, 0,'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "mso:f:p:qj:h", long_options, &option_index);
#else
        c = getopt(argc, argv, "mso:f:p:qj:h");
#endif
        if (c == -1)
            break;

        switch (c) {
        case 'p': // The precision (i.e. 32, 64, or arbitrary)
            if (std::string("32") == optarg) { precision = _4ti2_PREC_INT_32; }
            else if (std::string("64") == optarg) { precision = _4ti2_PREC_INT_64; }
            else if (std::string("arbitrary").find(optarg) == 0) { precision = _4ti2_PREC_INT_ARB; }
            else { 
                std::cerr << "Unrecognised option for -p,--precision: " << optarg << "\n";
                exit(1);
            }
            break;
        default:
            break;
        }
    }

    if (optind == argc-1) {
        filename = argv[optind];
    }
}
