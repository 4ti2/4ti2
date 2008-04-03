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

#include "../banner.h"

#include <iostream>
#include <fstream>

#include "Vector.hpp"
#include "VectorArray.hpp"
#include "Variables.hpp"
#include "Relation.hpp"
#include "LinearSystem.hpp"
#include "Lattice.hpp"
#include "Options.h"
#include "Exception.h"
#include "Algorithm.hpp"
#include "DefaultController.hpp"

using namespace _4ti2_zsolve_;

template <typename T> void read_relations (const Options& options, std::string name, LinearSystem <T>* system)
{
    std::ifstream file (name.c_str());
    if (!file.good ())
        return;

    size_t h,w;
    file >> h >> w;

    if (h != 1 || (w != system->relations () && system != NULL))
        throw IOException ("Badly formatted relations file " + name);

    if (system != NULL)
    {
        for (size_t i = 0; i < w; i++)
        {
            Relation <T> & rel = system->get_relation (i);
            char c;
            file >> c;
            switch (c)
            {
                case '<':
                    rel.set (Relation <T> :: LesserEqual);
                    if (options.graver () || options.hilbert ())
                        throw IOException ("Inequations for `graver` or `hilbert` executable. Use `zsolve` executable instead.");
                break;
                case '>':
                    rel.set (Relation <T> :: GreaterEqual);
                    if (options.graver () || options.hilbert ())
                        throw IOException ("Inequations for `graver` or `hilbert` executable. Use `zsolve` executable instead.");
                break;
                case '=':
                    rel.set (Relation <T> :: Equal);
                break;
                default:
                    std::string msg = "Unknown relation token '";
                    msg += c;
                    msg += "' in relations file ";
                    throw IOException (msg + name);
                break;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < w; i++)
        {
            char c;
            file >> c;
            if (c != '=')
                throw IOException ("Inequalities are not allowed, when reading from lattice.\n");
        }
    }
}

template <typename T> void read_signs (const Options& options, std::string name, VariableProperties <T>* properties)
{
    std::ifstream file (name.c_str());
    if (!file.good ())
        return;

    size_t h,w;
    file >> h >> w;
    if (h != 1 || w != properties->variables ())
        throw IOException ("Badly formatted sign file " + name);

    for (size_t i = 0; i < w; i++)
    {
        VariableProperty <T>& var = properties->get_variable (i);
        std::string s;
        file >> s;
        if (s == "free" || s == "f" || s == "0")
            var.set (true);
        else if (s == "hilbert" || s == "h" || s == "+h" || s == "1" || s == "+")
            var.set (false, 0, -1);
        else if (s == "-hilbert" || s == "-h" || s == "-1" || s == "-")
            var.set (false, 1, 0);
        else if (s == "graver" || s == "g" || s == "2")
        {
            if (options.hilbert ())
                throw IOException ("Graver components for `hilbert` executable. Use `graver` executable instead.");
            var.set (false);
        }
    }
}

template <typename T> void read_bounds (const Options& options, std::string name, VariableProperties <T>* properties, bool is_lower)
{
    std::ifstream file (name.c_str());
    if (!file.good ())
        return;

    size_t h, w;
    file >> h >> w;
    if (h != 1 || w != properties->variables ())
        throw IOException ("Badly formatted file " + name);

    for (size_t i = 0; i < w; i++)
    {
        VariableProperty <T>& var = properties->get_variable (i);
        std::string s;
        file >> s;
        if (s == "*")
        {
            if (is_lower && options.hilbert())
                throw IOException ("Lower bounds cannot be set to nonzero in `hilbert` executable. Use the `graver` executable instead.");
            var.set_bound (is_lower, is_lower ? 1 : -1);
        }
        else
        {
            T value = 0;
            parse_integer (s, value);
            if ((is_lower && value > 0) || (!is_lower && value < 0))
            {
                std::string first = (is_lower ? "Lower" : "Upper");
                first += " bound cannot be set to value ";
                std::string second = (is_lower ? "greater" : "lesser");
                second +=  " than zero.";
                throw IOException (first + second);
            }
            var.set_bound (is_lower, value);
        }
    }
}

template <typename T> Algorithm <T>* prepare_matrix (const Options& options, Controller <T> * controller, std::ifstream& matrix_file)
{
    // matrix
    VectorArray <T> matrix;
    matrix_file >> matrix;

    // rhs
    T* rhs;
    std::string rhs_name = options.project () + ".rhs";
    std::ifstream rhs_file (rhs_name.c_str ());
    if (rhs_file.good ())
    {
        size_t size;
        rhs_file >> size;
        if (size != 1)
            throw IOException ("Height of rhs must be 1!");
        rhs_file >> size;
        if (size != matrix.height ())
            throw IOException ("Height of matrix and size of right hand side differ!");
        rhs = read_vector <T> (rhs_file, size);
    }
    else
        rhs = create_zero_vector <T> (matrix.height ());

    if ((options.graver () || options.hilbert ()) && !is_zero_vector (rhs, matrix.height ()))
        throw IOException ("Right hand side must be a zero vector, when invoked as graver / hilbert binary!\n");

    LinearSystem <T> * system = new LinearSystem <T> (matrix, rhs, !(options.graver() || options.hilbert ()), options.hilbert () ? 0 : 1, -1);
    
    delete_vector (rhs);
    
    if (system->cancel_down ())
    {
        std::cout << "Canceled down linear system!\n" << std::endl;
    }

    read_relations <T> (options, options.project () + ".rel", system);
    read_signs <T> (options, options.project () + ".sign", system);
    read_bounds <T> (options, options.project () + ".lb", system, true);
    read_bounds <T> (options, options.project () + ".ub", system, false);

    Algorithm <T> * algorithm = new Algorithm <T> (system, controller);

    delete system;

    return algorithm;
}

template <typename T> Algorithm <T>* prepare_lattice (const Options& options, Controller <T> * controller, std::ifstream& lattice_file)
{
    // matrix
    VectorArray <T> vectors;
    lattice_file >> vectors;

    // rhs
    std::string rhs_name = options.project () + ".rhs";
    std::ifstream rhs_file (rhs_name.c_str ());
    if (rhs_file.good ())
    {
        T* rhs;
        size_t size;
        rhs_file >> size;
        if (size != 1)
            throw IOException ("Height of rhs must be 1!");
        rhs_file >> size;
        rhs = read_vector <T> (rhs_file, size);
        for (size_t i = 0; i < size; i++)
            if (rhs[i] != 0)
                throw IOException ("Right hand side must be a zero vector, when reading from lattice!\n");
        delete rhs;
    }

    bool free = false;
    T lower = 1;
    if (options.hilbert ())
        lower = 0;
    else if (!options.graver () && !options.hilbert ())
        free  = true;

    Lattice <T> * lattice = new Lattice <T> (&vectors, free, lower, -1);

    read_relations (options, options.project () + ".rel", (LinearSystem <T> *) NULL);
    read_signs (options, options.project () + ".sign", lattice);
    read_bounds (options, options.project () + ".lb", lattice, true);
    read_bounds (options, options.project () + ".ub", lattice, false);

    lattice->reduce_gaussian ();

    Algorithm <T> * algorithm = new Algorithm <T> (lattice, controller);

    delete lattice;

    return algorithm;
}

template <typename T> Algorithm <T>* prepare_backup (Options& options, Controller <T> * controller, std::ifstream& backup_file)
{
    backup_file >> options;

    Algorithm <T> * algorithm = new Algorithm <T> (backup_file, controller);

    return algorithm;
}

template <typename T> int zsolve_main (Options& options)
{
    Algorithm <T>* algorithm;
    DefaultController <T> * controller;
    std::ofstream* log_file = NULL;
    if (options.loglevel () > 0)
    {
        std::string log_name = options.project () + ".log";
        log_file = new std::ofstream (log_name.c_str(), options.resume () ? std::ios::out | std::ios::app : std::ios::out);
    }

    controller = new DefaultController <T> (&std::cout, log_file, options);

    std::string backup_name = options.project () + ".backup";
    std::ifstream backup_file (backup_name.c_str());

    if (backup_file.good ())
    {
        if (options.resume ())
            algorithm = prepare_backup <T> (options, controller, backup_file);
        else
            throw IOException ("Found backup file. Please restart with -r or delete backup file!\n", false);
    }
    else
    {
        if (options.resume ())
            throw IOException ("Started in resume mode, but no backup file found!\n", false);
        std::string matrix_name = options.project () + ".mat";
        std::ifstream matrix_file (matrix_name.c_str());
        
        if (matrix_file.good ())
            algorithm = prepare_matrix <T> (options, controller, matrix_file);
        else
        {
            std::string lattice_name = options.project () + ".lat";
            std::ifstream lattice_file (lattice_name.c_str());

            if (lattice_file.good ())
                algorithm = prepare_lattice <T> (options, controller, lattice_file);
            else
                throw IOException ("Neither " + options.project () + ".mat, " + options.project () + ".lat, nor " + options.project () + ".backup found!");
        }
    }

    algorithm->compute (options.backup_frequency ());

    algorithm->log_maxnorm ();

    if (options.graver ())
    {
        VectorArray <T> graver (algorithm->get_result_variables ());
        algorithm->extract_graver_results (graver);
        graver.save (options.project () + ".gra");
    }
    else if (options.hilbert ())
    {
        VectorArray <T> hilbert (algorithm->get_result_variables ());
        algorithm->extract_hilbert_results (hilbert);
        hilbert.save (options.project () + ".hil");
    }
    else
    {
        VectorArray <T> inhoms (algorithm->get_result_variables ());
        VectorArray <T> homs (algorithm->get_result_variables ());
        VectorArray <T> free (algorithm->get_result_variables ());

        algorithm->extract_zsolve_results (inhoms, homs, free);

        inhoms.save (options.project () + ".zinhom");
        homs.save (options.project () + ".zhom");
        if (free.vectors () > 0)
            free.save (options.project () + ".zfree");
    }

    delete algorithm;
    delete controller;

    if (log_file != NULL)
        delete log_file;

    return 0;
}

int main (int argc, char **argv)
{
    Options::print_banner ();

    Options options (argc, argv);
    if (options.verbosity () != 0)
    {
        options.print_precision ();
    }

    int result;
    try
    {
        if (options.precision () == 32)
            result = zsolve_main <int> (options);
        else if (options.precision () == 64)
            result = zsolve_main <long long> (options);
#ifdef _4ti2_GMP_
        else if (options.precision () == 0)
            result = zsolve_main <mpz_class> (options);
#endif
        else
            throw new IOException ("Invalid precision.");
    }
    catch (IOException e)
    {
        std::cerr << e << std::endl;
        if (e.print ())
            options.print_usage ();
        return 2;
    }

    return result;
}

