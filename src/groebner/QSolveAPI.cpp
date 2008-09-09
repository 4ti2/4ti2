#include <cstring>
#include <iostream>
#include <fstream>
#include "4ti2/4ti2.h"
#include "groebner/QSolveAPI.h"
#include "groebner/QSolveAlgorithm.h"
#include "groebner/VectorArrayAPI.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/LatticeBasis.h"

#include "Globals.h"
#include <unistd.h>

#ifdef _GNU_SOURCE
#include <getopt.h>
#endif

using namespace _4ti2_;

QSolveAPI::QSolveAPI()
{
    matrix = 0;
    sign = 0;
    rel = 0;
    lat = 0;
    ray = 0;
    cir = 0;
    qhom = 0;
    qfree = 0;
// TODO: We will need to change the next row.
#ifdef _4ti2_GMP_
    algorithm = SUPPORT;
#else
    algorithm = MATRIX;
#endif
    order = MAXCUTOFF;
}

QSolveAPI::~QSolveAPI()
{
    delete matrix;
    delete sign;
    delete rel;
    delete lat;
    delete ray;
    delete cir;
    delete qhom;
    delete qfree;
}

_4ti2_matrix*
QSolveAPI::create_matrix(int num_rows, int num_cols, const char* name)
{
    if (!strcmp(name, "mat")) { delete matrix; return (matrix = new VectorArrayAPI(num_rows, num_cols)); }
    if (!strcmp(name, "sign")) { delete sign; return (sign = new VectorArrayAPI(num_rows, num_cols)); }
    if (!strcmp(name, "rel")) { delete rel; return (rel = new VectorArrayAPI(num_rows, num_cols)); }
    if (!strcmp(name, "lat")) { delete lat; return (lat = new VectorArrayAPI(num_rows, num_cols)); }
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
    int m, n;
    in >> m >> n;
    _4ti2_matrix* matrix = create_matrix(m, n, name);
    matrix->read(in);
    return matrix;
}

_4ti2_matrix*
QSolveAPI::get_matrix(const char* name)
{
    if (!strcmp(name, "mat")) { return matrix; }
    if (!strcmp(name, "sign")) { return sign; }
    if (!strcmp(name, "rel")) { return rel; }
    if (!strcmp(name, "lat")) { return lat; }
    if (!strcmp(name, "ray")) { return ray; }
    if (!strcmp(name, "cir")) { return cir; }
    if (!strcmp(name, "qhom")) { return qhom; }
    if (!strcmp(name, "qfree")) { return qfree; }
    std::cerr << "ERROR: Unrecognised matrix type " << name << ".\n";
    return 0;
}

void
QSolveAPI::compute()
{
    // Consistency and default value checking.
    if (!matrix && !lat) {
        std::cerr << "ERROR: No matrix or lattice specified.\n";
        exit(1);
    }
    if (!matrix) {
        matrix = new VectorArrayAPI(0, lat->get_num_cols());
        lattice_basis(lat->data, matrix->data);
    }
    if (!lat) {
        lat = new VectorArrayAPI(0, matrix->get_num_cols());
        lattice_basis(matrix->data, lat->data);
    }
    if (!sign) {
        sign = new VectorArrayAPI(1, matrix->get_num_cols());
        for (Index i = 0; i < sign->get_num_cols(); ++i) { sign->data[0][i] = 0; }
    }
    if (!rel) {
        rel = new VectorArrayAPI(1, matrix->get_num_cols());
        for (Index i = 0; i < rel->get_num_cols(); ++i) { rel->data[0][i] = 0; }
    }
    assert(sign->get_number() == 1);
    assert(matrix->get_num_cols() == sign->get_num_cols());

    std::cout << "Matrix:\n";
    matrix->write(std::cout);
    std::cout << "Sign:\n";
    sign->write(std::cout);
    std::cout << "Rel:\n";
    rel->write(std::cout);
    std::cout << "Lat:\n";
    lat->write(std::cout);

    // Delete previous computation.
    delete ray; delete cir; delete qhom; delete qfree;
    ray = new VectorArrayAPI(0, matrix->get_num_cols());
    ray->data.insert(lat->data);
    cir = new VectorArrayAPI(0, matrix->get_num_cols());
    qhom = new VectorArrayAPI(0, matrix->get_num_cols());
    qfree = new VectorArrayAPI(0, matrix->get_num_cols());

    QSolveAlgorithm alg(algorithm, order);
    alg.compute(matrix->data, ray->data, cir->data, qfree->data, rel->data[0], sign->data[0]); 

    ray->data.sort();
    cir->data.sort();
    qfree->data.sort();
    VectorArray::transfer(ray->data, qhom->data);
    VectorArray cir_data_neg(cir->data);
    VectorArray::transfer(cir->data, qhom->data);
    cir_data_neg.mul(-1);
    VectorArray::transfer(cir_data_neg, qhom->data);
}

void
QSolveAPI::set_options(int argc, char** argv)
{
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
            {"quiet",        0, 0,'q'},
            {"help",         0, 0,'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "mso:f:p:qh",
                 long_options, &option_index);
#else
        c = getopt(argc, argv, "mso:f:p:qh");
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

    if (optind != argc) {
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
  PROJECT.mat         A matrix (optional if lattice basis is given).\n\
  PROJECT.lat         A lattice basis (optional if matrix is given).\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative, '0' means a free variable, and '2' means\n\
                      both non-negative and non-positive).\n\
                      It is optional, and the default is all free.\n\
  PROJECT.rel         The relations on the matrix rows ('<','>','=').\n\
                      It is optional and the default is all '='.\n\
                      The matrix must be given with this file.\n";
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
  -m, --matrix               Use the Matrix algorithm (default for 32 and 64).\n\
  -s, --support              Use the Support algorithm (default for arbitrary).\n\
  -o, --order=ORDERING       Set ORDERING as the ordering in which the columns\n\
                             are chosen. The possible orderings are `maxinter',\n\
                             `minindex', `maxcutoff' (default), and `mincutoff'.\n\
  -f, --output_freq=n        Set the frequency of output (default is 1000).\n\
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
    // First, we get rid of any previous data structures.
    delete matrix; delete sign; delete rel; delete lat;
    matrix = 0; sign = 0; rel = 0; lat = 0;

    std::string basename(basename_c_str);   

    // Read in the file with the matrix.
    std::string matrix_filename(basename + ".mat");
    create_matrix(matrix_filename.c_str(), "mat");
    if (matrix == 0) {
        create_matrix(basename.c_str(), "mat");
        if (matrix != 0) {
            std::cout << "WARNING: Please specify the matrix in '" << matrix_filename;
            std::cout << "' instead of '" << basename << "'.\n";
        }
    }

    // Read in the file with the sign constraints.
    std::string sign_filename(basename + ".sign");
    create_matrix(sign_filename.c_str(), "sign");

    // Read in the file with the matrix relations.
    std::string rel_filename(basename + ".rel");
    create_matrix(rel_filename.c_str(), "rel");

    // Read in the file with the matrix relations.
    std::string ray_filename(basename + ".lat");
    create_matrix(ray_filename.c_str(), "lat");
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
