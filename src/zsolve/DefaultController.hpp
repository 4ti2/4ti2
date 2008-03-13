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

#ifndef _4ti2_zsolve__DefaultController_
#define _4ti2_zsolve__DefaultController_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "zsolve/Options.h"
#include "zsolve/Controller.hpp"
#include "zsolve/Timer.h"
#include "zsolve/Algorithm.hpp"

template <typename T> class DefaultController : public Controller <T>
{
protected:
    std::ostream* m_console;
    std::ofstream* m_log;
    const Options& m_options;
    Timer m_all_timer;
    Timer m_var_timer;
    Timer m_sum_timer;
    Timer m_norm_timer;
    T m_last_maxnorm;
    size_t m_last_maxnorm_vectors;

public:
    DefaultController (std::ostream* console, std::ofstream* log, const Options& options) : m_options (options)
    {
        m_console = console;
        m_log = log;
	m_last_maxnorm = 0;
	m_last_maxnorm_vectors = 0;
    }

    ~DefaultController () {}

    void log_system (LinearSystem <T> * system)
    {
        if (m_options.verbosity () != 0)
            *m_console << "Linear system to solve:\n\n" << *system << std::endl;
        if (m_options.loglevel () > 0)
            *m_log << "Linear system to solve:\n\n" << *system << std::endl;
    }

    void log_homogenized_system (LinearSystem <T> * system)
    {
        if (m_options.verbosity () != 0)
            *m_console << "Linear system of homogeneous equalities to solve:\n\n" << *system << std::endl;
        if (m_options.loglevel () > 0)
            *m_log << "Linear system of homogeneous equalities to solve:\n\n" << *system << std::endl;
    }

    void log_lattice (Lattice <T> * lattice)
    {
        if (m_options.verbosity () != 0)
            *m_console << "Lattice:\n\n" << *lattice << std::endl;
        if (m_options.loglevel () > 0)
            *m_log << "Lattice:\n\n" << *lattice << std::endl;
    }

    void log_variable_start (size_t variable, size_t vectors)
    {
        m_var_timer.reset();

        if (m_options.verbosity () == 1)
        {
            *m_console << "Appending variable " << variable << " ..." << std::flush;
        }
        else if (m_options.verbosity () > 1)
        {
            if (variable > 1)
                *m_console << '\n';
            *m_console << "Appending variable " << variable << ".\n" <<std::endl;
        }

        if (m_options.loglevel () == 1)
        {
            *m_log << "Appending variable " << variable << " ..." << std::flush;
        }
        else if (m_options.loglevel () > 1)
        {
            if (variable > 1)
                *m_log << '\n';
            *m_log << "Appending variable " << variable << ".\n" <<std::endl;
        }
    }

    void log_variable_end (size_t variable, size_t vectors)
    {
        if (m_options.verbosity () == 1)
        {
            *m_console << " Solutions: " << vectors << ", Step: " << m_var_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        }
        else if (m_options.verbosity () > 1)
        {
            if (m_options.verbosity () == 2)
                *m_console << "\n";
            *m_console << "Finished variable " << variable << ". Solutions: " << vectors << ", Step: " << m_var_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        }

        if (m_options.loglevel () == 1)
        {
            *m_log << " Solutions: " << vectors << ", Step: " << m_var_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        }
        else if (m_options.loglevel () > 1)
        {
            if (m_options.loglevel () == 2)
                *m_log << "\n";
            *m_log << "Finished variable " << variable << ". Solutions: " << vectors << ", Step: " << m_var_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        }
    }

    void log_sum_start (size_t variable, const T& sum, size_t vectors)
    {
        m_sum_timer.reset ();

        if (m_options.verbosity () == 2)
            *m_console << "  Variable: " << variable << ", Sum = " << sum << " ..." << std::flush;
        else if (m_options.verbosity () == 3)
            *m_console << "  Variable: " << variable << ", Processing sum " << sum << "\n" << std::endl;
        if (m_options.loglevel () == 2)
            *m_log << "  Variable: " << variable << ", Sum = " << sum << " ..." << std::flush;
        else if (m_options.loglevel () == 3)
            *m_log << "  Variable: " << variable << ", Processing sum " << sum << "\n" << std::endl;
    }

    void log_sum_end (size_t variable, const T& sum, size_t vectors) 
    {
        if (m_options.verbosity () == 2)
            *m_console << " Solutions: " << vectors << ", Step: " << m_sum_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        else if (m_options.verbosity () == 3)
            *m_console << "\n  Finished sum " << sum << ". Solutions: " << vectors << ", Step: " << m_sum_timer << "s, Time: " << m_all_timer << "s\n" << std::endl;
        if (m_options.loglevel () == 2)
            *m_log << " Solutions: " << vectors <<  ", Step: " << m_sum_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        else if (m_options.loglevel () == 3)
            *m_log << "\n  Finished sum " << sum << ". Solutions: " << vectors << ", Step: " << m_sum_timer << "s, Time: " << m_all_timer << "s\n" << std::endl;
    }

    void log_norm_start (size_t variable, const T& sum, const T& norm, size_t vectors)
    {
        m_norm_timer.reset ();

        if (m_options.verbosity () == 3)
            *m_console << "    Variable: " << variable << ", Norm = " << norm << " + " << sum-norm << " ..." << std::flush;
        if (m_options.loglevel () == 3)
            *m_log << "    Variable: " << variable << ", Norm = " << norm << " + " << sum-norm << " ..." << std::flush;
    }

    void log_norm_end (size_t variable, const T& sum, const T& norm, size_t vectors)
    {
        if (m_options.verbosity () == 3)
            *m_console << " Solutions: " << vectors << ", Step: " << m_norm_timer << "s, Time: " << m_all_timer << "s" << std::endl;
        if (m_options.loglevel () == 3)
            *m_log << " Solutions: " << vectors << ", Step: " << m_norm_timer << "s, Time: " << m_all_timer << "s" << std::endl;
    }

    void log_result (size_t inhoms, size_t homs, size_t free)
    {
        if (m_options.verbosity () != 0)
            *m_console << "\nFinal basis has " << inhoms << " inhomogeneous, " << homs << " homogeneous and " << free << " free elements. Time: " << m_all_timer << "s" << std::endl;
        if (m_options.loglevel () != 0)
            *m_log << "\nFinal basis has " << inhoms << " inhomogeneous, " << homs << " homogeneous and " << free << " free elements. Time: " << m_all_timer << "s" << std::endl;
    }

    void log_status (size_t variable, const T& sum, const T& max_sum, const T& norm, size_t vectors, int backup_frequency, Timer& timer)
    {
        static int wrap = 1000;
        if (m_options.verbosity () < 0)
        {
            static int i = 0;
            if (i == 0)
                i = wrap;
            i--;
            if (i==1)
            {
                static Timer wrap_timer;
                // update status every 1-2 seconds
                if (wrap_timer.get_elapsed_time () > 1.0)
                    wrap /= 2;
                else if (wrap_timer.get_elapsed_time () < 0.5)
                    wrap *= 2;

                // get needed space
                static unsigned int max_space = 0;
                std::stringstream ss;

                if (m_options.verbosity () == -1)
                {
                    ss << "\rVariable: " << variable << ", Sum: " << sum << ", Norm: " << norm << ", Max-Norm: " << max_sum << ", Solutions: " << vectors << ", Time: " << m_all_timer << "s" << std::flush;
                }
                else if (backup_frequency != 0)
                {
                    double next_backup = backup_frequency - timer.get_elapsed_time ();
                    ss << "\rVariable: " << variable << ", Sum: " << sum << ", Norm: " << norm << " + " << (sum-norm) << ", Solutions: " << vectors;
                    ss << ", Time (norm): " << m_norm_timer << "s, Time (sum): " << m_sum_timer << "s, Time (variable): " << m_var_timer << "s, Time: " << m_all_timer << "s, Next backup: ";
                    if (next_backup >= 0.0)
                        ss << next_backup << "s" << std::flush;
                    else
                        ss << "on next step" << std::flush;
                }
                else
                {
                    ss << "\rVariable: " << variable << ", Sum: " << sum << ", Norm: " << norm << " + " << (sum-norm) << ", Max-Norm: " << max_sum << ", Solutions: " << vectors;
                    ss << ", Time (norm): " << m_norm_timer << "s, Time (sum): " << m_sum_timer << "s, Time (variable): " << m_var_timer << "s, Time: " << m_all_timer << "s" << std::flush;
                }
                
                // fill with whitespace
                std::string s = ss.str ();
                std::string fill;
                if (s.size () > max_space)
                    max_space = s.size ();
                else
                    for (unsigned int i = s.size (); i < max_space; i++)
                        fill = fill + " ";
                *m_console << s << fill << std::flush;
                *m_console << s << std::flush;
                wrap_timer.reset ();
            }
        }
    }
    
    void log_resume (size_t variables, size_t variable, const T& sum, const T& norm, size_t vectors)
    {
        if (m_options.verbosity () != 0)
            *m_console << "Resuming backup at variable " << variable << " of " << variables << ", sum " << sum << " (" << norm << " + " << (sum - norm) << ")" << ", with " << vectors << " solutions.\n" << std::endl;
        if (m_options.loglevel () != 0)
            *m_log << "\n\nResuming backup at variable " << variable << " of " << variables << ", sum " << sum << " (" << norm << " + " << (sum - norm) << ")" << ", with " << vectors << " solutions.\n" << std::endl;
    }

    void log_maxnorm (Algorithm <T> * algorithm, bool final)
    {
        if (m_options.maxnorm () == 1 && final || m_options.maxnorm () == 2)
        {
            VectorArray <T> maxnorm (algorithm->get_result_variables ());
            T norm = algorithm->extract_maxnorm_results (maxnorm);

            if (final)
            {
                if (m_options.verbosity () != 0)
                    *m_console << "\nFinal basis has " << maxnorm.vectors () << " vectors with a maximum norm of " << norm << "." << std::endl;
                if (m_options.loglevel () != 0)
                    *m_log << "\nFinal basis has " << maxnorm.vectors () << " vectors with a maximum norm of " << norm << "." << std::endl;
                maxnorm.save (m_options.project () + ".maxnorm");
            }
            else if (norm != m_last_maxnorm || maxnorm.vectors () != m_last_maxnorm_vectors)
            {
                m_last_maxnorm_vectors = maxnorm.vectors ();
                m_last_maxnorm = norm;
                if (m_options.verbosity () != 0)
                    *m_console << "\nBasis has " << m_last_maxnorm_vectors << " vectors with a maximum norm of " << norm << "." << std::endl;
                if (m_options.loglevel () != 0)
                    *m_log << "\nBasis has " << m_last_maxnorm_vectors << " vectors with a maximum norm of " << norm << "." << std::endl;
                maxnorm.save (m_options.project () + ".maxnorm");
            }
        }
    }
    
    void save_lattice (Lattice <T> * lattice)
    {
        std::string name = m_options.project () + ".lat";
        std::ofstream file (name.c_str(), std::ios::out);

        VectorArray <T> * va = lattice;
        file << *va << std::endl;
    }

    void backup_data (Lattice <T> & lattice, size_t current, const T& sum, const T& norm, bool symmetric)
    {
        // write to backup~ - otherwise, everything can get lost, if one does abort the program while backing up.
        std::string name = m_options.project () + ".backup~";
        std::ofstream file (name.c_str(), std::ios::out);

        // OPTIONS

        // verbosity, loglevel, backup frequency
        file << m_options.verbosity () << "\n";
        file << m_options.loglevel () << "\n";
        file << m_options.backup_frequency () << "\n";
        // graver, hilbert, zsolve
        if (m_options.graver ())
            file << "g\n";
        else if (m_options.hilbert ())
            file << "h\n";
        else
            file << "z\n";
        // maxnorm
        file << (m_options.maxnorm () ? "1\n" : "0\n");
        // precision
        if (m_options.precision () == 32)
            file << "32\n";
        else if (m_options.precision () == 64)
            file << "64\n";
        else
            file << "gmp\n";

        // CONTROLLER DATA
        file << "\n";

        // times

        file << m_all_timer.get_elapsed_time () << " " << m_var_timer.get_elapsed_time () << " " << m_sum_timer.get_elapsed_time () << "\n";

        // ALGORITHM DATA
        file << "\n";

        // current, sum, norm, symmetric
        file << current << " " << sum << " " << norm << " " << (symmetric ? "1 " : "0 ") << "\n";
        file << (int) lattice.vectors () << " " << (int) lattice.variables () << "\n";
        for (size_t i = 0; i < lattice.variables (); i++)
        {
            const VariableProperty <T>& var =  lattice.get_variable (i);
            var.dump (file);
            file << "\n";
        }
        for (size_t i = 0; i < lattice.vectors (); i++)
        {
            print_vector (file, lattice[i], lattice.variables ());
            file << "\n";
        }
        file << std::flush;
        file.close ();

        // now move the backup~ to backup
        std::string new_name = m_options.project () + ".backup";
        rename (name.c_str(), new_name.c_str());

        if (m_options.verbosity () > 0)
           *m_console << " Paused for backup.\nResuming computation ..." << std::flush;
        if (m_options.loglevel () > 0)
           *m_log << " Paused for backup.\nResuming computation ..." << std::flush;
    }

    void read_backup (std::ifstream& in)
    {
        in >> m_all_timer >> m_var_timer >> m_sum_timer;
    }
};

#endif
