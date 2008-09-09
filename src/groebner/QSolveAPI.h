#ifndef _4ti2_groebner__QSolveAPI_
#define _4ti2_groebner__QSolveAPI_

#include "groebner/4ti2API.h"
#include "groebner/QSolveVariant.h"
#include "groebner/QSolveConsOrder.h"

namespace _4ti2_ {

class VectorArrayAPI;

class QSolveAPI : public _4ti2_state {
public:
    QSolveAPI();
    virtual ~QSolveAPI();

    virtual void compute();

    virtual void set_options(int argc, char** argv);

    virtual void read(const char* basename);
    virtual void write(const char* basename);

    virtual _4ti2_matrix* create_matrix(int num_rows, int num_cols, const char* name);
    virtual _4ti2_matrix* create_matrix(const char* filename, const char* name);
    virtual _4ti2_matrix* create_matrix(std::istream& in, const char* name);

    virtual _4ti2_matrix* get_matrix(const char* name);

protected:
    QSolveVariant algorithm;
    QSolveConsOrder order;

    virtual void write_usage();
    virtual void write_options();
    virtual void write_input_files();
    virtual void write_output_files();

    void unrecognised_option_argument(const char* option);

    VectorArrayAPI* matrix;
    VectorArrayAPI* sign;
    VectorArrayAPI* rel;
    VectorArrayAPI* lat;
    VectorArrayAPI* ray;
    VectorArrayAPI* cir;
    VectorArrayAPI* qhom;
    VectorArrayAPI* qfree;
};

} // namespace _4ti2_

#endif
