/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2008 4ti2 team.
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

#ifndef _4ti2_zsolve__ZSolveAPI_
#define _4ti2_zsolve__ZSolveAPI_

#include <cstring>

#include "4ti2/4ti2xx.h"
#include "zsolve/VectorArrayAPI.hpp"
#include "zsolve/RelAPI.h"
#include "zsolve/SignAPI.h"
#include "zsolve/BoundAPI.hpp"

#include "zsolve/Options.h"
#include "zsolve/Vector.hpp"
#include "zsolve/VectorArray.hpp"
#include "zsolve/Variables.hpp"
#include "zsolve/Relation.hpp"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Lattice.hpp"
#include "zsolve/Options.h"
#include "zsolve/Exception.h"
#include "zsolve/Algorithm.hpp"
#include "zsolve/DefaultController.hpp"

namespace _4ti2_zsolve_ {

template <class T>
class ZSolveAPI : public _4ti2_state {
public:
    ZSolveAPI();
    virtual ~ZSolveAPI();

    virtual void compute();

    virtual void set_options(int argc, char** argv);
    virtual void set_options(const Options& o);

    virtual void read(const char* project);
    virtual void write(const char* project);

    virtual _4ti2_matrix* create_matrix(int num_rows, int num_cols, const char* name);
    virtual _4ti2_matrix* create_matrix(const char* filename, const char* name);
    virtual _4ti2_matrix* create_matrix(std::istream& in, const char* name);

    virtual _4ti2_matrix* get_matrix(const char* name);

protected:

    // Extract the output after running the algorithm.
    virtual void extract_results(Algorithm <T>* algorithm);

    // Checks whether the input data is in a consistent state.
    virtual void check_consistency();

    Options options;

    VectorArrayAPI<T>* mat;
    VectorArrayAPI<T>* lat;
    VectorArrayAPI<T>* rhs;
    BoundAPI<T>* ub;
    BoundAPI<T>* lb;
    RelAPI* rel;
    SignAPI* sign;
    VectorArrayAPI<T>* zinhom;
    VectorArrayAPI<T>* zhom;
    VectorArrayAPI<T>* zfree;

    bool free_default;
    T lower_default;
    T upper_default;
};

template <class T>
ZSolveAPI<T>::ZSolveAPI()
        : mat(0), lat(0), rhs(0), ub(0), lb(0), rel(0), sign(0), zinhom(0), zhom(0), zfree(0),
          free_default(true), lower_default(1), upper_default(-1)
{
}

template <class T>
ZSolveAPI<T>::~ZSolveAPI()
{
    delete mat; mat = 0;
    delete lat; lat = 0;
    delete rhs; rhs = 0;
    delete ub; ub = 0;
    delete lb; lb = 0;
    delete rel; rel = 0;
    delete sign; sign = 0;
    delete zinhom; zinhom = 0;
    delete zhom; zhom = 0;
    delete zfree; zfree = 0;
}

template <class T>
void
ZSolveAPI<T>::set_options(int argc, char** argv)
{
    options.process_options(argc, argv);
}

template <class T>
void
ZSolveAPI<T>::set_options(const Options& o)
{
    options = o;
}

template <class T>
void
ZSolveAPI<T>::read(const char* project_c_str)
{
    std::string project(project_c_str);   

    // Read in all the possible input files.
    create_matrix((project + ".mat").c_str(), "mat");
    create_matrix((project + ".lat").c_str(), "lat");
    create_matrix((project + ".rhs").c_str(), "rhs");
    create_matrix((project + ".ub").c_str(), "ub");
    create_matrix((project + ".lb").c_str(), "lb");
    create_matrix((project + ".rel").c_str(), "rel");
    create_matrix((project + ".sign").c_str(), "sign");
}

template <class T>
void
ZSolveAPI<T>::check_consistency()
{
    if (!mat && !lat) {
        throw IOException("No `mat' or `lat' specified!");
    }
    if (mat && lat) {
        throw IOException("Both `mat' and `lat' cannot be given as input!");
    }

    if (lat && rhs) {
        throw IOException ("Both `lat' and `rhs' cannot be given as input!");
    }
    if (lat && rel) {
        throw IOException ("Both `lat' and `rel' cannot be given as input!");
    }

    int num_variables = 0;
    if (mat) { num_variables = mat->get_num_cols(); }
    if (lat) { num_variables = lat->get_num_cols(); }

    if (rhs && rhs->get_num_rows() != 1) {
        throw IOException ("Height of `rhs' should be 1!");
    }
    if (mat && rel && rel->get_num_cols() != mat->get_num_rows()) {
        throw IOException ("Height of `mat' and size of `rel' differ!");
    }
    if (mat && rhs && rhs->get_num_cols() != mat->get_num_rows()) {
        throw IOException ("Height of `mat' and size of `rhs' differ!");
    }
    if (ub && ub->get_num_cols() != num_variables) {
        throw IOException ("Width of `mat' and size of `ub' differ!");
    }
    if (lb && lb->get_num_cols() != num_variables) {
        throw IOException ("Width of `mat' and size of `lb' differ!");
    }
    if (sign && sign->get_num_cols() != num_variables) {
        throw IOException ("Width of `mat' and size of `sign' differ!");
    }
}

template <class T>
void
ZSolveAPI<T>::write(const char* project_c_str)
{
    std::string project(project_c_str);   

    if (zinhom) { zinhom->write((project + ".zinhom").c_str()); }
    if (zhom) { zhom->write((project + ".zhom").c_str()); }
    if (zfree && zfree->data.height() > 0) { zfree->write((project + ".zfree").c_str()); }
}

template <class T>
_4ti2_matrix*
ZSolveAPI<T>::create_matrix(int num_rows, int num_cols, const char* name)
{
    if (!strcmp(name, "mat")) { delete mat; return (mat = new VectorArrayAPI<T>(num_rows, num_cols)); }
    if (!strcmp(name, "lat")) { delete lat; return (lat = new VectorArrayAPI<T>(num_rows, num_cols)); }
    if (!strcmp(name, "rhs")) { delete rhs; return (rhs = new VectorArrayAPI<T>(num_rows, num_cols)); }
    if (!strcmp(name, "lb")) { delete lb; return (lb = new BoundAPI<T>(num_rows, num_cols, true)); }
    if (!strcmp(name, "ub")) { delete ub; return (ub = new BoundAPI<T>(num_rows, num_cols, false)); }
    if (!strcmp(name, "rel")) { delete rel; return (rel = new RelAPI(num_rows, num_cols)); }
    if (!strcmp(name, "sign")) { delete sign; return (sign = new SignAPI(num_rows, num_cols)); }
    std::cerr << "ERROR: Unrecognised input matrix type " << name << ".\n";
    return 0;
}

template <class T>
_4ti2_matrix*
ZSolveAPI<T>::create_matrix(const char* filename, const char* name)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    return create_matrix(file, name);
}

template <class T>
_4ti2_matrix*
ZSolveAPI<T>::create_matrix(std::istream&in, const char* name)
{
    int m, n;
    in >> m >> n;
    _4ti2_matrix* mat = create_matrix(m, n, name);
    mat->read(in);
    return mat;
}

template <class T>
_4ti2_matrix*
ZSolveAPI<T>::get_matrix(const char* name)
{
    if (!strcmp(name, "mat")) { return mat; }
    if (!strcmp(name, "lat")) { return lat; }
    if (!strcmp(name, "rhs")) { return rhs; }
    if (!strcmp(name, "ub")) { return ub; }
    if (!strcmp(name, "lb")) { return lb; }
    if (!strcmp(name, "rel")) { return rel; }
    if (!strcmp(name, "sign")) { return sign; }
    if (!strcmp(name, "zhom")) { return zhom; }
    if (!strcmp(name, "zinhom")) { return zinhom; }
    if (!strcmp(name, "zfree")) { return zfree; }
    std::cerr << "ERROR: Unrecognised mat type " << name << ".\n";
    return 0;
}

template <class T>
void
ZSolveAPI<T>::compute()
{
    check_consistency();
    
    Algorithm <T>* algorithm;
    DefaultController <T> * controller;
    std::ofstream* log_file = 0;
    if (options.loglevel () > 0) {
        std::string log_name = options.project () + ".log";
        log_file = new std::ofstream (log_name.c_str(), options.resume () ? std::ios::out | std::ios::app : std::ios::out);
    }
    controller = new DefaultController <T> (&std::cout, log_file, options);

    std::string backup_name = options.project () + ".backup";
    std::ifstream backup_file (backup_name.c_str());

    if (backup_file.good () && !options.resume ()) {
        throw IOException ("Found backup file. Please restart with -r or delete backup file!\n", false);
    }

    if (options.resume()) {
        if (!backup_file.good()) {
            throw IOException ("Started in resume mode, but no backup file found!\n", false);
        }
        backup_file >> options;
        algorithm = new Algorithm <T> (backup_file, controller);
    }
    else if (mat) {
        // TODO: transfer rhs, ub, lb, sign and rel.
        T* rhs_vec = create_zero_vector <T> (mat->data.height());
        if (rhs) { 
            for (size_t i = 0; i < rhs->data.width(); ++i) {
                rhs_vec[i] = rhs->data[0][i];
            }
        }
        LinearSystem <T> * system = new LinearSystem <T> (mat->data, rhs_vec,
                        free_default, lower_default, upper_default);
        delete_vector(rhs_vec);
        if (sign) {
            for (size_t i = 0; i < sign->data.width(); ++i) {
                switch (sign->data[0][i]) {
                    case 0:
                        system->get_variable(i).set(true);
                        break;
                    case 1:
                        system->get_variable(i).set(false, 0, -1);
                        break;
                    case -1:
                        system->get_variable(i).set(false, 1, 0);
                        break;
                    case 2:
                        system->get_variable(i).set(false);
                        break;
                    default:
                        // TODO: The following error message should be more informative.
                        throw IOException("Unknown sign value.");
                }
            }
        }
        if (rel) {
            for (size_t i = 0; i < rel->data.width(); ++i) {
                switch (rel->data[0][i]) {
                    case 0:
                        system->get_relation(i).set(Relation<T> :: Equal);
                        break;
                    case 1:
                        system->get_relation(i).set(Relation<T> :: GreaterEqual);
                        break;
                    case -1:
                        system->get_relation(i).set(Relation<T> :: LesserEqual);
                        break;
                    default:
                        // TODO: The following error message should be more informative.
                        throw IOException("Unknown relation value.");
                }
            }
        }
        if (lb) {
            for (size_t i = 0; i < lb->data.width(); ++i) {
                system->get_variable(i).set_bound(true, lb->data[0][i]);
            }
        }
        if (ub) {
            for (size_t i = 0; i < ub->data.width(); ++i) {
                system->get_variable(i).set_bound(false, ub->data[0][i]);
            }
        }

        system->cancel_down();
        algorithm = new Algorithm <T> (system, controller);
        delete system;
    }
    else if (lat) {
        // TODO: transfer ub, lb, and sign.
        Lattice <T> * lattice = new Lattice <T> (&lat->data, free_default, lower_default, upper_default);
        if (sign) {
            for (size_t i = 0; i < sign->data.width(); ++i) {
                switch (sign->data[0][i]) {
                    case 0:
                        lattice->get_variable(i).set(true);
                        break;
                    case 1:
                        lattice->get_variable(i).set(false, 0, -1);
                        break;
                    case -1:
                        lattice->get_variable(i).set(false, 1, 0);
                        break;
                    case 2:
                        lattice->get_variable(i).set(false);
                        break;
                    default:
                        // TODO: The following error message should be more informative.
                        throw IOException("Unknown sign value.");
                }
            }
        }
        if (lb) {
            for (size_t i = 0; i < lb->data.width(); ++i) {
                lattice->get_variable(i).set_bound(true, lb->data[0][i]);
            }
        }
        if (ub) {
            for (size_t i = 0; i < ub->data.width(); ++i) {
                lattice->get_variable(i).set_bound(false, ub->data[0][i]);
            }
        }

        lattice->reduce_gaussian();
        algorithm = new Algorithm <T> (lattice, controller);
        delete lattice;
    }
    else {
        throw IOException ("Neither " + options.project () + ".mat, " + options.project () + ".lat, nor "
                        + options.project () + ".backup found!");
    }

	algorithm->compute (options.backup_frequency ());

    algorithm->log_maxnorm ();

    extract_results(algorithm);

    delete algorithm;
    delete controller;

    if (log_file) { delete log_file; }
}

template <class T>
void
ZSolveAPI<T>::extract_results(Algorithm <T>* algorithm)
{
    delete zinhom;
    delete zhom;
    delete zfree;
    zinhom = new VectorArrayAPI <T> (0, algorithm->get_result_variables ());
    zhom = new VectorArrayAPI <T> (0, algorithm->get_result_variables ());
    zfree = new VectorArrayAPI <T> (0, algorithm->get_result_variables ());
    algorithm->extract_zsolve_results (zinhom->data, zhom->data, zfree->data);
}

} // namespace _4ti2_

#endif
