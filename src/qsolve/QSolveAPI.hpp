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

#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "qsolve/QSolveAPI.h"
#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/MatrixAlgorithm.h"
#include "qsolve/SupportAlgorithm.h"
#include "qsolve/VectorArray.h"
#include "qsolve/VectorArrayAPI.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/Cone.h"
#include "qsolve/IndexSetD.h"
#include "qsolve/IndexSetDS.h"

#include "qsolve/Globals.h"
#include "qsolve/Debug.h"

#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

using namespace _4ti2_;

template <class T>
QSolveAPI<T>::QSolveAPI()
{
    sign_default = _4ti2_FR;
    rel_default = _4ti2_EQ;

#ifdef _4ti2_GMP_
    algorithm = SUPPORT;
#else
    algorithm = MATRIX;
#endif
    order = MAXCUTOFF;
}

template <class T>
QSolveAPI<T>::~QSolveAPI()
{
}

template <class T>
_4ti2_matrix*
QSolveAPI<T>::create_matrix(int m, int n, const char* name)
{
    if (!strcmp(name, "mat")) { mat.init(m,n); return &mat; }
    if (!strcmp(name, "sign")) { sign.init(m,n); return &sign; }
    if (!strcmp(name, "rel")) { rel.init(m,n); return &rel; }
    std::cerr << "ERROR: Unrecognised input matrix type " << name << ".\n";
    return 0;
}

template <class T>
_4ti2_matrix*
QSolveAPI<T>::create_matrix(const char* filename, const char* name)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    return create_matrix(file, name);
}

template <class T>
_4ti2_matrix*
QSolveAPI<T>::create_matrix(std::istream&in, const char* name)
{
    int m, n;
    m = n = 0;
    in >> m >> n;
    // TODO: Input error detection.
    _4ti2_matrix* mat = create_matrix(m, n, name);
    mat->read(in);
    return mat;
}

template <class T>
_4ti2_matrix*
QSolveAPI<T>::get_matrix(const char* name)
{
    if (!strcmp(name, "mat")) { return &mat; }
    if (!strcmp(name, "sign")) { return &sign; }
    if (!strcmp(name, "rel")) { return &rel; }
    if (!strcmp(name, "ray")) { return &ray; }
    if (!strcmp(name, "cir")) { return &cir; }
    if (!strcmp(name, "qhom")) { return &qhom; }
    if (!strcmp(name, "qfree")) { return &qfree; }
    std::cerr << "ERROR: Unrecognised mat type " << name << ".\n";
    return 0;
}

template <class T>
void
QSolveAPI<T>::pre_compute()
{
    Size n = mat.get_num_cols();
    Size m = mat.get_num_rows();

    if (sign.get_num_rows() != 0) {
        if (sign.get_num_rows() != 1 || sign.get_num_cols() != n) { // TODO: Use exceptions.
            std::cerr << "ERROR: Incorrect sizes for sign matrix.\n";
            exit(1);
        }
        for (Index i = 0; i < n; ++i) {
            mat.cone.set_constraint_type(i, (_4ti2_constraint) sign.data[0][i]);
        }
    }
    else {
        for (Index i = 0; i < n; ++i) { 
            mat.cone.set_constraint_type(i, sign_default);
        }
    }

    if (rel.get_num_rows() != 0) {
        if (rel.get_num_rows() != 1 || rel.get_num_cols() != m) { // TODO: Use exceptions.
            std::cerr << "ERROR: Incorrect sizes for rel matrix.\n";
            exit(1);
        }
        for (Index i = 0; i < m; ++i) {
            mat.cone.set_constraint_type(i+n, (_4ti2_constraint) rel.data[0][i]);
        }
    }
    else {
        for (Index i = n; i < n+m; ++i) {
            mat.cone.set_constraint_type(i, rel_default);
        }
    }

    // Delete previous computation.
    ray.init(0, mat.get_num_cols());
    cir.init(0, mat.get_num_cols());
    qhom.init(0, mat.get_num_cols());
    qfree.init(0, mat.get_num_cols());
}

template <class T>
void
QSolveAPI<T>::compute()
{
    print_banner();
    pre_compute();

    Size n = mat.get_num_cols();
    Size m = mat.get_num_rows();
    // TODO: Change the following size for the IndexSet.
    if (n+m <= IndexSetDS::max_size) {
        QSolveAlgorithm<T,IndexSetDS>* alg = 0;
        if (algorithm == MATRIX) { alg = new MatrixAlgorithm<T,IndexSetDS>(order); }
        else { alg = new SupportAlgorithm<T,IndexSetDS>(order); }
        alg->compute(mat.cone, ray.data, cir.data, qfree.data);
        delete alg;
    } else {
        QSolveAlgorithm<T,IndexSetD>* alg = 0;
        if (algorithm == MATRIX) { alg = new MatrixAlgorithm<T,IndexSetD>(order); }
        else { alg = new SupportAlgorithm<T,IndexSetD>(order); }
        alg->compute(mat.cone, ray.data, cir.data, qfree.data);
        delete alg;
    }

    post_compute();
}


template <class T>
void
QSolveAPI<T>::post_compute()
{
    ray.data.sort();
    qfree.data.sort();
    qhom.data.transfer(ray.data, 0, ray.data.get_number(), qhom.data.get_number());
    VectorArrayT<T> cir_data_neg(cir.data);
    cir_data_neg.mul(-1);
    qhom.data.transfer(cir.data, 0, cir.data.get_number(), qhom.data.get_number());
    qhom.data.transfer(cir_data_neg, 0, cir_data_neg.get_number(), qhom.data.get_number());
    qhom.data.sort();
}

template <class T>
void
QSolveAPI<T>::set_options(int argc, char** argv)
{
    int c;
    optind = 1;
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

        c = getopt_long (argc, argv, "mso:f:p:qj:h",
                 long_options, &option_index);
#else
        c = getopt(argc, argv, "mso:f:p:qj:h");
#endif
        if (c == -1)
            break;

        switch (c) {
        case 'm':
            algorithm = MATRIX;
            break;
        case 's':
            algorithm = SUPPORT;
            break;
        case 'o':
            if (std::string("maxinter").find(optarg) == 0)
            { order = MAXINTER; }
            else if (std::string("minindex").find(optarg) == 0)
            { order = MININDEX; }
            else if (std::string("maxcutoff").find(optarg) == 0)
            { order = MAXCUTOFF; }
            else if (std::string("mincutoff").find(optarg) == 0)
            { order = MINCUTOFF; }
            else { unrecognised_option_argument("-o, --order"); }
            break;
        case 'q':
            out = new std::ofstream;
            //TODO: err = new std::ofstream;
            break;
        case 'f':
            if (sscanf(optarg, "%d", &Globals::output_freq) != 1)
            {  unrecognised_option_argument("-f, --output_freq"); }
            break;
        case 'p': // The precision (i.e. 32, 64, or arbitrary)
            if (std::string("32").find(optarg) == 0) { }
            else if (std::string("64").find(optarg) == 0) { }
            else if (std::string("arbitrary").find(optarg) == 0) { }
            else { unrecognised_option_argument("-p, --precision"); }
            break;
        case 'j': // Multithreaded run.
            if (sscanf(optarg, "%d", &Globals::num_threads) != 1)
            {  unrecognised_option_argument("-j, --threads"); }
            break;
        case 'h':
        case '?':
        case ':':
            write_usage();
            exit(1);
            break;

        default:
            std::cerr << "ERROR: getopt returned unknown character code";
            std::cerr << std::endl;
            write_usage();
            exit(1);
        }
    }

    if (optind != argc && optind != argc-1) {
        std::cerr << "ERROR: unrecognised options ... ";
        for (; optind < argc; ++optind) {
            std::cerr << " " << argv[optind];
        }
        std::cerr << "\n";
        write_usage();
        exit(1);
    }
}

template <class T>
void
QSolveAPI<T>::write_usage()
{
    std::cerr << "Usage: qsolve [options] <PROJECT>\n\n";
    std::cerr << "Computes a generator description of a cone.\n";
    write_input_files();
    write_output_files();
    write_options();
}

template <class T>
void
QSolveAPI<T>::write_input_files()
{
    std::cerr << "\
Input Files:\n\
  PROJECT.mat         A matrix (compulsory).\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative, '0' means a free variable, and '2' means\n\
                      both non-negative and non-positive).\n\
                      It is optional, and the default is all free.\n\
  PROJECT.rel         The relations on the matrix rows ('<','>','=').\n\
                      It is optional and the default is all '='.\n\
                      The mat must be given with this file.\n";
}

template <class T>
void
QSolveAPI<T>::write_output_files()
{
    std::cerr << "\
Output Files:\n\
  PROJECT.qhom        The homogeneous generators of the linear system.\n\
  PROJECT.qfree       A basis for the linear subspace of the cone.\n\
                      If this file does not exist then the linear subspace \n\
                      is trivial.\n\n";
}

template <class T>
void
QSolveAPI<T>::write_options()
{
    std::cerr << "\
Options:\n\
  -p, --precision=PREC       Select PREC as the integer arithmetic precision.\n\
                             PREC is one of the following: `64' (default),\n\
                             `32', and `arbitrary' (only `arb` is needed).\n\
  -m, --mat               Use the Matrix algorithm (default for 32 and 64).\n\
  -s, --support              Use the Support algorithm (default for arbitrary).\n\
  -o, --order=ORDERING       Set ORDERING as the ordering in which the columns\n\
                             are chosen. The possible orderings are `maxinter',\n\
                             `minindex', `maxcutoff' (default), and `mincutoff'.\n\
  -f, --output_freq=n        Set the frequency of output (default is 1000).\n\
  -j, --threads=n            Set the maximum number of threads to n.\n\
  -q, --quiet                Do not output anything to the screen.\n\
  -h, --help                 Display this help and exit.\n\
\n\
Only short options are supported on sun machines.\n\
\n";
}

template <class T>
void
QSolveAPI<T>::unrecognised_option_argument(const char* option)
{
   std::cerr << "4ti2: ";
   std::cerr << "Unrecognised argument \"" << optarg << "\" ";
   std::cerr << "for the option " << option << ".\n\n";
   write_usage();
   exit(1);
}

template <class T>
void
QSolveAPI<T>::read(const char* basename_c_str)
{
    std::string basename(basename_c_str);   

    // Read in the file with the mat.
    std::string mat_filename(basename + ".mat");
    if (create_matrix(mat_filename.c_str(), "mat") == 0) {
        if (create_matrix(basename.c_str(), "mat") == 0) {
            std::cerr << "ERROR: No constraint matrix specified.\n";
            std::cerr << "ERROR: Expected matrix in '" << mat_filename << "'\n";
            exit(1);
        }
        else {
            std::cerr << "WARNING: Please specify the matrix in '" << mat_filename;
            std::cerr << "' instead of '" << basename << "'.\n";
        }
    }

    // Read in the file with the sign constraints.
    std::string sign_filename(basename + ".sign");
    create_matrix(sign_filename.c_str(), "sign");

    // Read in the file with the mat relations.
    std::string rel_filename(basename + ".rel");
    create_matrix(rel_filename.c_str(), "rel");
}

template <class T>
void
QSolveAPI<T>::write(const char* basename_c_str)
{
    std::string basename(basename_c_str);

    std::string qhom_filename(basename + ".qhom");
    qhom.write(qhom_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    qfree.write(qfree_filename.c_str());
}

template <class T>
void 
QSolveAPI<T>::print_banner()
{
    *out << FORTY_TWO_BANNER;
}
