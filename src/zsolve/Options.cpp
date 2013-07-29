/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter.

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

#include <cstdlib>
#include <cstring>
#include "../banner.h"
#include "Options.h"

#ifdef __GNU_LIBRARY__
#	include <getopt.h>
#endif

namespace _4ti2_zsolve_ {

Options::Options ()
{
    set_defaults ();
}

Options::Options (const Options& o)
{
    *this = o;
}

Options::Options (int argc, char **argv)
{
    process_options(argc, argv);
}

Options&
Options::operator= (const Options& o)
{
    m_project = o.m_project;
    m_graver = o.m_graver;
    m_hilbert = o.m_hilbert;
    m_precision = o.m_precision;
    m_verbosity = o.m_verbosity;
    m_loglevel = o.m_loglevel;
    m_backup_frequency = o.m_backup_frequency;
    m_resume = o.m_resume;
    m_maxnorm = o.m_maxnorm;
    return *this;
}

void
Options::set_defaults ()
{
    m_project = "zsolve";
    m_graver = false;
    m_hilbert = false;
    m_precision = _4ti2_PREC_INT_32;
    m_verbosity = -1;
    m_loglevel = 0;
    m_backup_frequency = 0;
    m_resume = false;
    m_maxnorm = false;
}

void
Options::process_options (int argc, char** argv)
{
    set_defaults();
    optind = 1;

    int c;

#ifdef __GNU_LIBRARY__
	static struct option long_options[] =
	{
		{ "backup", optional_argument, NULL, 'b'},
		{ "help", no_argument, NULL, 'h'},
		{ "log", optional_argument, NULL, 'l'},
		{ "quiet", no_argument, NULL, 'q'},
		{ "resume", no_argument, NULL, 'r'},
		{ "update", optional_argument, NULL, 'u'},
		{ "verbose", optional_argument, NULL, 'v'},
		{ "hilbert", no_argument, NULL, 'H'},
		{ "graver", no_argument, NULL, 'G'},
		{ "maxnorm", no_argument, NULL, 'm'},
		{ "precision", required_argument, NULL, 'p'}
	};
	while ((c = getopt_long(argc, argv, "b::hl::qru::v::HGmp:", long_options, NULL)) != -1)
#else
	while ((c = getopt(argc, argv, "b::hl::qru::v::HGmp:")) != -1)
#endif
	{
		if (optarg != NULL && optarg[0] == '=')
			optarg++;
		switch (c)
		{
			case 'b':
				if (optarg == NULL)
                    m_backup_frequency = 3600;
				else
				{
                    char d;
                    int b;
					if (sscanf (optarg, "%d%c", &b, &d) == 2)
					{
						if (d == 'm')
                            m_backup_frequency = b * 60;
						else if (d == 'h')
							m_backup_frequency = b * 3600;
						else if (d == 'd')
							m_backup_frequency = b * 86400;
						else if (d != 's')
						{
                            std::cout << "Invalid backup modifier: " << (char) d << std::endl;
                            exit (1);
						}
                        else
                            m_backup_frequency = b;
					}
					else
					{
                        std::cout << "Invalid backup argument: -b" << optarg << std::endl;
						exit(1);
					}
				}
			break;
            case 'h':
                print_usage ();
                exit (0);
            break;
			case 'l':
				if (optarg == NULL || !strcmp (optarg, "1"))
                    m_loglevel = 1;
				else if (!strcmp (optarg, "2") || !strcmp (optarg, "l"))
					m_loglevel = 2;
				else if (!strcmp (optarg, "3") || !strcmp (optarg, "ll"))
					m_loglevel = 3;
				else if (!strcmp (optarg, "0"))
					m_loglevel = 0;
				else
				{
                    std::cout << "Invalid loglevel: -v" << optarg << std::endl;
					exit(1);
				}
            break;
            case 'q':
                m_verbosity = 0;
			break;
            case 'r':
                m_resume = true;
            break;
            case 'u':
				if (optarg == NULL || !strcmp (optarg, "1"))
                    m_verbosity = -1;
				else if (!strcmp (optarg, "2") || !strcmp (optarg, "u"))
                    m_verbosity = -2;
				else if (!strcmp (optarg, "0"))
                    m_verbosity = 0;
				else
				{
                    std::cout << "Invalid verbosity/update: -u" << optarg << std::endl;
					exit(1);
				}
            break;
            case 'v':
				if (optarg == NULL || !strcmp (optarg, "1"))
                    m_verbosity = 1;
				else if (!strcmp (optarg, "2") || !strcmp (optarg, "v"))
                    m_verbosity = 2;
				else if (!strcmp (optarg, "3") || !strcmp (optarg, "vv"))
                    m_verbosity = 3;
				else if (!strcmp (optarg, "0"))
                    m_verbosity = 0;
				else
				{
                    std::cout << "Invalid verbosity/update: -v" << optarg << std::endl;
					exit(1);
				}
            break;
			case 'H':
				m_hilbert = true;
			break;
			case 'G':
                m_graver = true;
			break;
			case 'm':
				m_maxnorm = true;
			break;
			case 'p':
                if (optarg == NULL || !strcmp (optarg, "32"))
                    m_precision = _4ti2_PREC_INT_32;
                else if (!strcmp (optarg, "64"))
                    m_precision = _4ti2_PREC_INT_64;
                else if (!strcmp (optarg, "gmp") || !strcmp (optarg, "arbitrary"))
                {
                    m_precision = _4ti2_PREC_INT_ARB;
#ifndef _4ti2_HAVE_GMP
                    std::cout << "This binary was compiled without GMP support!" << std::endl;
                    exit (1);
#endif
                }
                else
                {
                    std::cout << "Invalid precision: -p" << optarg << std::endl;
                    exit (1);
                }
			break;
			default:
                std::cout << "Unknown getopt sequence " << c << ", " << optarg << std::endl;
                print_usage ();
                exit (1);
			break;
		}
	}

    if (m_hilbert && m_graver)
	{
        std::cout << "Input error: A combination of -H and -G is not allowed!" << std::endl;
		exit(1);
	}

    if (optind == argc-1) { m_project = argv[optind]; }
	else if (optind > argc) { print_usage (); exit (1); }
	else if (optind < argc-1) 
    {
        std::cerr << "Argument error: Only one project file is possible: You specified '";
        std::cerr << argv[optind] << "' and '" << argv[optind+1] << "'!\n";
		exit(1);
	}
}

void Options::print_banner ()
{
    std::cout << FORTY_TWO_BANNER << "\n" << std::flush;
}
    
void Options::print_usage () const
{
    std::cout << "Usage: ";
    if (m_graver)
        std::cout << "graver";
    else if (m_hilbert)
        std::cout << "hilbert";
    else
        std::cout << "zsolve";
    std::cout << " [options] PROJECT\n\n";

    std::cout << "[Basic options]\n\n";

#ifdef _4ti2_HAVE_GMP
    std::cout << " -p=PREC, --precision=PREC  Use precision (32, 64, gmp). Default is 32 bit\n";
#else
    std::cout << " -p=PREC, --precision=PREC  Use precision (32, 64). Default is 32 bit\n";
#endif
    std::cout << " -m, --maxnorm              Write vectors with maximum norm to PROJECT.maxnorm\n";
    std::cout << " -b[FREQ], --backup[=FREQ]  Frequently backup status to PROJECT.backup\n";
    std::cout << " -r, --resume               Resume from backup file PROJECT.backup\n";
    std::cout << " -h, --help                 Display this help\n";
    std::cout << "\n";

    std::cout << "[Output options]\n\n";
    std::cout << " -q, --quiet        Quit mode\n";
    std::cout << " -u, --update[=1]   Updated output on console (default)\n";
    std::cout << " -uu, --update=2    More verbose updated output on console\n";
    std::cout << " -v, --verbose[=1]  Output once every variable computation\n";
    std::cout << " -vv, --verbose=2   Output once every norm sum computation\n";
    std::cout << " -vvv, --verbose=3  Output once every norm computation\n";
    std::cout << "\n";

    std::cout << "[Logging options]\n\n";
    std::cout << " -n, --log=0    Disable logging (default)\n";
    std::cout << " -l, --log[=1]  Log once every variable computation to PROJECT.log\n";
    std::cout << " -ll, --log=2   Log once every norm sum computation to PROJECT.log\n";
    std::cout << " -lll, --log=3  Log once every norm computation to PROJECT.log\n";
    std::cout << "\n";

    std::cout << "[Used files]\n\n";
    std::cout << "PROJECT.mat     Matrix\n";
    std::cout << "PROJECT.rhs     Right hand side\n";
    std::cout << "PROJECT.rel     Relations (<, >, =)\n";
    std::cout << "PROJECT.sign    Sign of columns (optional)\n";
    std::cout << "PROJECT.lb      Lower bounds of columns (optional)\n";
    std::cout << "PROJECT.ub      Upper bounds of columns (optional)\n";
    std::cout << "\n";
    std::cout << "PROJECT.backup  Backup file\n";
    std::cout << "PROJECT.backup~ Temporary backup file - if it exsts, it may be newer than PROJECT.backup!\n";
    std::cout << "\n";
    std::cout << "PROJECT.zinhom  Inhomogeneous part of the solution\n";
    std::cout << "PROJECT.zhom    Homogeneous part of the solution\n";
    std::cout << "PROJECT.zfree   Free part of the solution\n";
    std::cout << "PROJECT.maxnorm Vectors with maximum norm\n";
    std::cout << std::endl;
}
    
void Options::print_precision () const
{
    if (m_precision == _4ti2_PREC_INT_32)
        std::cout << "Using 32 bit integers.\n" << std::endl;
    else if (m_precision == _4ti2_PREC_INT_64)
        std::cout << "Using 64 bit integers.\n" << std::endl;
    else
        std::cout << "Using arbitrary precision integers.\n" << std::endl;
}

std::string Options::project () const
{
    return m_project;
}

int Options::verbosity () const
{
    return m_verbosity;
}

int Options::loglevel () const
{
    return m_loglevel;
}

int Options::backup_frequency () const
{
    return m_backup_frequency;
}

bool Options::resume () const
{
    return m_resume;
}

bool Options::hilbert () const
{
    return m_hilbert;
}

bool Options::graver () const
{
    return m_graver;
}

bool Options::maxnorm () const
{
    return m_maxnorm;
}

_4ti2_precision Options::precision () const
{
    return m_precision;
}

std::istream& operator>>(std::istream& in, Options& options)
{
    int v,l,m;
    int b;
    std::string mode;
    std::string prec;

    in >> v >> l >> b >> mode >> m >> prec;

    if (v != options.m_verbosity)
    {
        std::cout << "Option warning: Verbosity from backup file, line 1 (" << v << ") and command line option (" << options.verbosity () << ") differ!\n" << std::endl;
    }
    if (l != options.m_loglevel)
    {
        std::cout << "Option warning: Loglevel from backup file, line 2 (" << l << ") and command line option (" << options.loglevel () << ") differ!\n" << std::endl;
    }
    if (options.m_backup_frequency == 0)
    {
        std::cout << "Option error: Backup is deactivated for resume. If you really like to do this, please change line 3 of " << options.project () << ".backup!\n" << std::endl;
        exit (1);
    }
    if ((mode == "g" && !options.m_graver) || (mode == "h" && !options.m_hilbert) || (mode == "z" && (options.m_graver || options.m_hilbert)))
    {
        std::cout << "Option error: Mode (graver, hilbert, zsolve) from backup file, line 4 (" << mode << ") and command line option differ!\n If you like to change it for resume, edit the backup file!\n" << std::endl;
        exit (1);
        
    }
    if ((prec == "32" && options.m_precision != _4ti2_PREC_INT_32) || (prec == "64" && options.m_precision != _4ti2_PREC_INT_64) || (prec == "gmp" && options.m_precision == _4ti2_PREC_INT_ARB))
    {
        std::cout << "Option error: Precision from backup file, line 6 (" << prec << ") and command line option (";
        if (options.precision () == 0)
            std::cout << "gmp";
        else
            std::cout << options.precision ();
        std::cout <<  ") differ!\n If you like to change it for resume, edit the backup file!\n" << std::endl;
        exit (1);
    }

    return in;
}

} // namespace _4ti2_zsolve_
