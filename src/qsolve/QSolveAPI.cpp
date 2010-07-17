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
#include <sstream>
#include <unistd.h>
#include <algorithm>

#include "qsolve/QSolveAPI.h"
#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/MatrixAlgorithm.h"
#include "qsolve/SupportAlgorithm.h"
#include "qsolve/VectorArray.h"
#include "qsolve/MatrixAPI.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/Cone.h"
#include "qsolve/IndexSetD.h"
#include "qsolve/IndexSetDS.h"
#include "qsolve/MatrixAPI.h"
#include "qsolve/MatrixWrapper.h"

#include "qsolve/Globals.h"
#include "qsolve/Debug.h"

#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

using namespace _4ti2_;

void
QSolveAPI::initialise_data()
{
    switch (prec) {
    case _4ti2_PREC_INT_ARB:
#ifdef _4ti2_GMP_
        initialise_dataT<mpz_class>(); break;
#else
        std::cerr << "ERROR: Arbitrary precision not supported.\n";        
        exit(1);
#endif
    case _4ti2_PREC_INT_32:
        initialise_dataT<int32_t>(); break;
    case _4ti2_PREC_INT_64:
        initialise_dataT<int64_t>(); break;
    }

    if (prec == _4ti2_PREC_INT_ARB) { algorithm = SUPPORT; }
    else { algorithm = MATRIX; }
}

template <class T>
void
QSolveAPI::initialise_dataT()
{
    delete mat; delete ray; delete cir; delete qhom; delete qfree;
    mat = new MatrixWrapper<VectorArrayT<T> >();
    ray = new MatrixAPI<VectorArrayT<T> >();
    cir = new MatrixAPI<VectorArrayT<T> >();
    qhom = new MatrixAPI<VectorArrayT<T> >();
    qfree = new MatrixAPI<VectorArrayT<T> >();
}

QSolveAPI::QSolveAPI()
    : _4ti2_state(), prec(_4ti2_PREC_INT_64), input(DEFAULT), mat(0), ray(0), cir(0), qhom(0), qfree(0)
{
    print_banner();
    sign_default = _4ti2_FR;
    rel_default = _4ti2_EQ;

    sign = new MatrixAPI<VectorArrayT<int32_t> >;
    rel = new MatrixAPI<VectorArrayT<int32_t> >;
    initialise_data();
    order = MAXCUTOFF;
}

QSolveAPI::~QSolveAPI()
{
    delete mat; delete sign; delete rel; delete ray; delete cir; delete qhom; delete qfree;
}

_4ti2_matrix*
QSolveAPI::create_matrix(int m, int n, const char* name)
{
    if (!strcmp(name, "mat")) { mat->resize(m,n); return mat; }
    if (!strcmp(name, "sign")) { sign->resize(m,n); return sign; }
    if (!strcmp(name, "rel")) { rel->resize(m,n); return rel; }
    std::cerr << "ERROR: Unrecognised input matrix type " << name << ".\n";
    return 0;
}

_4ti2_matrix*
QSolveAPI::create_matrix(const char* filename, const char* name)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    return create_matrix(file, name);
}

_4ti2_matrix*
QSolveAPI::create_matrix(std::istream&in, const char* name)
{
    int m,n;
    m = n = 0;
    bool readline = true;
    while (readline) {
        in >> m >> n;
        if (in.good()) { readline = false; } // Stop if ok.
        if (in.eof()) { return 0; }
        in.clear();
        while (in.get() != '\n') { if (!in.good()) { return 0; } }
    }
    _4ti2_matrix* matrix = create_matrix(m, n, name);
    matrix->read(in);
    return matrix;
}

_4ti2_matrix*
QSolveAPI::get_matrix(const char* name)
{
    if (!strcmp(name, "mat")) { return mat; }
    if (!strcmp(name, "sign")) { return sign; }
    if (!strcmp(name, "rel")) { return rel; }
    if (!strcmp(name, "ray")) { return ray; }
    if (!strcmp(name, "cir")) { return cir; }
    if (!strcmp(name, "qhom")) { return qhom; }
    if (!strcmp(name, "qfree")) { return qfree; }
    std::cerr << "ERROR: Unrecognised mat type " << name << ".\n";
    return 0;
}

void
QSolveAPI::pre_compute()
{
    Size n = mat->get_num_cols();
    Size m = mat->get_num_rows();
    if (sign->get_num_rows() != 0) {
        if (sign->get_num_rows() != 1 || sign->get_num_cols() != n) { // TODO: Use exceptions.
            std::cerr << "ERROR: Incorrect sizes for sign matrix.\n";
            exit(1);
        }
    }
    else {
        sign->resize(1,n);
        for (Index i = 0; i < n; ++i) { 
            sign->set_entry(0, i, (int32_t) sign_default);
        }
    }

    if (rel->get_num_rows() != 0) {
        if (rel->get_num_rows() != 1 || rel->get_num_cols() != m) { // TODO: Use exceptions.
            std::cerr << "ERROR: Incorrect sizes for rel matrix.\n";
            exit(1);
        }
    }
    else {
        rel->resize(1,m);
        for (Index i = 0; i < m; ++i) {
            rel->set_entry(0, i, (int32_t) rel_default);
        }
    }

    // Delete previous computation.
    ray->resize(0, mat->get_num_cols());
    cir->resize(0, mat->get_num_cols());
    qhom->resize(0, mat->get_num_cols());
    qfree->resize(0, mat->get_num_cols());
}

void
QSolveAPI::compute()
{
    pre_compute();

    switch (prec) {
    case _4ti2_PREC_INT_ARB:
        *out << "Using Abitrary precision.\n";
#ifdef _4ti2_GMP_
        computeT<mpz_class>(); break;
#else
        std::cerr << "ERROR: Arbitrary precision not supported.\n";        
        exit(1);
#endif
    case _4ti2_PREC_INT_32:
        *out << "Using 32 bit precision.\n";
        computeT<int32_t>(); break;
    case _4ti2_PREC_INT_64:
        *out << "Using 64 bit precision.\n";
        computeT<int64_t>(); break;
    }

    post_compute();
}

template <class T>
void
QSolveAPI::computeT()
{
    Size n = mat->get_num_cols();
    Size m = mat->get_num_rows();
    ConeT<T> cone(dynamic_cast<MatrixWrapper<VectorArrayT<T> >&>(*mat));
    int32_t c;
    for (Index i = 0; i < n; ++i) { 
        sign->get_entry(0,i,c);
        cone.set_constraint_type(i,(_4ti2_constraint) c);
    }
    for (Index i = n; i < n+m; ++i) { 
        rel->get_entry(0,i-n,c);
        cone.set_constraint_type(i,(_4ti2_constraint) c);
    }
    VectorArrayT<T>& ray_data = dynamic_cast<MatrixAPI<VectorArrayT<T> >&>(*ray).data;
    VectorArrayT<T>& cir_data = dynamic_cast<MatrixAPI<VectorArrayT<T> >&>(*cir).data;
    VectorArrayT<T>& qfree_data = dynamic_cast<MatrixAPI<VectorArrayT<T> >&>(*qfree).data;

    QSolveAlgorithm<T>* alg = 0;
    if (algorithm == MATRIX) { alg = new MatrixAlgorithm<T>(order); }
    else { alg = new SupportAlgorithm<T>(order); }
    alg->compute(cone, ray_data, cir_data, qfree_data);
    delete alg;
}

void
QSolveAPI::post_compute()
{
    //ray->data.sort();
    //qfree->data.sort();
    qhom->swap(*ray);
#if 0
    VectorArrayT<T> cir_data_neg(*cir);
    cir_data_neg.mul(-1);
    qhom->transfer(*cir, 0, cir->get_num_rows(), qhom->get_num_rows());
    qhom->transfer(cir_data_neg, 0, cir_data_neg.get_number(), qhom->data.get_number());
    qhom->data.sort();
#endif
}

void
QSolveAPI::set_options(int argc, char** argv)
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
            {"input",        1, 0,'i'},
            {"threads",      1, 0,'j'},
            {"newline",      0, 0,'n'},
            {"quiet",        0, 0,'q'},
            {"help",         0, 0,'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "mso:f:p:i:nqj:h",
                 long_options, &option_index);
#else
        c = getopt(argc, argv, "mso:f:p:i:nqj:h");
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
            if (std::string("32").find(optarg) == 0) { prec = _4ti2_PREC_INT_32; }
            else if (std::string("64").find(optarg) == 0) { prec = _4ti2_PREC_INT_64; }
            else if (std::string("arbitrary").find(optarg) == 0) { prec = _4ti2_PREC_INT_ARB; }
            else { unrecognised_option_argument("-p, --precision"); }
            initialise_data();
            break;
        case 'i':
            if (std::string("ine") == optarg) { input = INE; }
            else if (std::string("ext") == optarg) { input = EXT; }
            else if (std::string("ieq") == optarg) { input = IEQ; }
            else if (std::string("poi") == optarg) { input = POI; }
            else { unrecognised_option_argument("-i, --input"); }
        case 'j': // Multithreaded run.
            if (sscanf(optarg, "%d", &Globals::num_threads) != 1)
            {  unrecognised_option_argument("-j, --threads"); }
            break;
        case 'n':
            ENDL = '\n';
            break;
        case 'h':
        case '?':
        case ':':
            write_usage();
            exit(1);
            break;

        default:
            std::cerr << "ERROR: getopt returned unknown character code.\n";
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

void
QSolveAPI::write_usage()
{
    std::cerr << "Usage: qsolve [options] <PROJECT>\n\n";
    std::cerr << "Computes a generator description of a cone.\n";
    write_input_files();
    write_output_files();
    write_options();
}

void
QSolveAPI::write_input_files()
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

void
QSolveAPI::write_output_files()
{
    std::cerr << "\
Output Files:\n\
  PROJECT.qhom        The homogeneous generators of the linear system.\n\
  PROJECT.qfree       A basis for the linear subspace of the cone.\n\
                      If this file does not exist then the linear subspace \n\
                      is trivial.\n\n";
}

void
QSolveAPI::write_options()
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
  -j, --threads=n            Set the maximum number of threads to n. Default is 1.\n\
  -q, --quiet                Do not output anything to the screen.\n\
  -h, --help                 Display this help and exit.\n\
\n\
Only short options are supported on sun machines.\n\
\n";
}

void
QSolveAPI::unrecognised_option_argument(const char* option)
{
   std::cerr << "4ti2: ";
   std::cerr << "Unrecognised argument \"" << optarg << "\" ";
   std::cerr << "for the option " << option << ".\n\n";
   write_usage();
   exit(1);
}

void
QSolveAPI::read(const char* basename_c_str)
{
    std::string basename(basename_c_str);   

    // Check to see whether there is a cdd ine input file.
    if (basename.size() > 4 && basename.compare(basename.size()-4,4,".ine") == 0) {
        *out << "Using cdd ine input file " << basename << ".\n";
        if (!create_matrix(basename.c_str(), "mat")) {
            std::cerr << "Could not parse file " << basename << ".\n";
            return;
        }
        int n = mat->get_num_cols();
        int m = mat->get_num_rows();
        sign->resize(1,n);
        for (Index i = 0; i < n; ++i) { 
            sign->set_entry(0, i, (int32_t) _4ti2_FR);
        }
        // The variable corresponding to the rhs is >= 0.
        sign->set_entry(0, 0, (int32_t) _4ti2_LB);
        rel->resize(1,m);
        for (Index i = 0; i < m; ++i) { 
            rel->set_entry(0, i, (int32_t) _4ti2_LB);
        }
        return;
    }

    // Check to see whether there is a cdd ext input file.
    if (basename.size() > 4 && basename.compare(basename.size()-4,4,".ext") == 0) {
        *out << "Using cdd ext input file " << basename << ".\n";
        if (!create_matrix(basename.c_str(), "mat")) {
            std::cerr << "Could not parse file " << basename << ".\n";
            return;
        }
        int n = mat->get_num_cols();
        int m = mat->get_num_rows();
        sign->resize(1,n);
        for (Index i = 0; i < n; ++i) { 
            sign->set_entry(0, i, (int32_t) _4ti2_FR);
        }
        rel->resize(1,m);
        for (Index i = 0; i < m; ++i) { 
            rel->set_entry(0, i, (int32_t) _4ti2_LB);
        }
        return;
    }

    // Check to see whether there is a porta ieq input file.
    if (basename.size() > 4 && basename.compare(basename.size()-4,4,".ieq") == 0) {
        std::ifstream file(basename.c_str());
        parse_porta_ieq(file);
        return;
    }

    // Read in the file with the matrix.
    std::string mat_filename(basename + ".mat");
    if (create_matrix(mat_filename.c_str(), "mat") != 0 
        || create_matrix(basename.c_str(), "mat") != 0) {
        // Read in the file with the sign constraints.
        std::string sign_filename(basename + ".sign");
        create_matrix(sign_filename.c_str(), "sign");

        // Read in the file with the mat relations.
        std::string rel_filename(basename + ".rel");
        create_matrix(rel_filename.c_str(), "rel");
        return;
    }

    std::cerr << "ERROR: No constraint matrix specified.\n";
    exit(1); // TODO
}

void
QSolveAPI::write(const char* basename_c_str)
{
    std::string basename(basename_c_str);

    std::string qhom_filename(basename + ".qhom");
    qhom->write(qhom_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    qfree->write(qfree_filename.c_str());
}

#if 0
void
QSolveAPI::parse_project_file(std::istream& in)
{
    int m,n;
    m = n = 0;
    std::string s;
    while (readline) {
    
    }
}

_4ti2_matrix*
QSolveAPI::create_matrix(std::istream&in, const char* name)
{
    int m,n;
    m = n = 0;
    bool readline = true;
    while (readline) {
        in >> m >> n;
        if (in.good()) { readline = false; } // Stop if ok.
        if (in.eof()) { return 0; }
        in.clear();
        while (in.get() != '\n') { if (!in.good()) { return 0; } }
    }
    _4ti2_matrix* matrix = create_matrix(m, n, name);
    matrix->read(in);
    return matrix;
}
#endif

void
QSolveAPI::parse_porta_ieq(std::istream& in)
{
    std::string input;
    getline(in, input);
    // strip whitespace
    input.erase(std::remove_if(input.begin(), input.end(), isspace), input.end());
    std::cout << input << std::endl;
    int n = -1;
    sscanf(input.c_str(), "DIM=%d", &n);
    if (n == -1) { 
        std::cerr << "INPUT ERROR: did not find DIM in " << basename << "\n";
        exit(1); // TODO: Use exceptions.
    }
    ++n; // Allow for an extra entry for the rhs.
    while (in.good()) {
        in.clear();
        getline(in, input);
        // strip whitespace
        input.erase(std::remove_if(input.begin(), input.end(), isspace), input.end());
        std::cout << input << std::endl;
        if (input == "INEQUALITIES_SECTION") { break; }
    }
    if (!in.good()) { 
        std::cerr << "INPUT ERROR: did not find INEQUALITIES_SECTION in " << basename << "\n";
        exit(1); // TODO: Use exceptions.
    }

    int m = 0;
    create_matrix(m, n, "mat");
    std::vector<_4ti2_constraint> cons;
    input.clear();
    while (!getline(in, input).eof()) {
        // strip input of whitespace.
        input.erase(std::remove_if(input.begin(), input.end(), isspace), input.end());
        if (input == "END") { break; }
        if (input.empty()) { continue; } // skip empty lines.
        mat->resize(m+1,n); // Add an extra row.
        for (int i = 0; i < n; ++i) { mat->set_entry(m, i, 0); } // Zero row.
        std::istringstream ss(input);
        char c;
        int coeff, index, sign;
        while (ss.good()) {
            sign = coeff = 1;
            ss >> c;
            if (c == '<' || c == '>' || c == '=') { break; }
            if (c == '+') { sign = 1; ss >> c; } else if (c == '-') { sign = -1; ss >> c; }
            if (isdigit(c)) {
                ss.unget();
                ss >> coeff >> c;
            }
            if (!isalpha(c)) { std::cerr << "INPUT ERROR: " << ss.str() << "\n"; exit(1); }
            coeff *= sign;
            ss >> index;
            if (!ss.good()) { std::cerr << "INPUT ERROR: " << ss.str() << "\n"; exit(1); }
            std::cout << std::showpos << coeff << std::noshowpos << c << index << " ";
            mat->set_entry(m, index, coeff);
        }
        if (!ss.good()) { std::cerr << "INPUT ERROR: " << ss.str() << "\n"; exit(1); }
        std::cout << c;
        _4ti2_constraint type = _4ti2_EQ;
        if (c == '=') { type = _4ti2_EQ; ss >> c; }
        if (c == '<') { type = _4ti2_UB; ss >> c; }
        else if (c == '>') { type = _4ti2_LB; ss >> c; }
        if (isdigit(c)) { ss.unget(); }
        else { std::cout << c; }
        ss >> coeff;
        if (ss.fail()) { std::cerr << "INPUT ERROR: " << ss.str() << "\n"; exit(1); }
        std::cout << " " << coeff << std::endl;
        std::cout << input << std::endl;
        mat->set_entry(m, 0, coeff);
        cons.push_back(type);
        ++m;
        input.clear();
    }

    rel->resize(1,m);
    for (int i = 0; i < m; ++i) {
        rel->set_entry(0, i, (int32_t) cons[i]);
    }
    sign->resize(1,n);
    for (int i = 0; i < n; ++i) {
        sign->set_entry(0, i, (int32_t) _4ti2_FR);
    }
    mat->write(std::cout);
    rel->write(std::cout);
}

void 
QSolveAPI::print_banner()
{
    *out << FORTY_TWO_BANNER;
}
