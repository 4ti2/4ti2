#include "4ti2/4ti2.h"
#include "groebner/4ti2API.h"
#include "groebner/QSolveAPI.h"
#include "groebner/RaysAPI.h"
#include "groebner/CircuitsAPI.h"

// TODO: Handle different precision.
// TODO: Handle errors.

using namespace _4ti2_;

extern "C" 
{

// Global variable (Argh!) representing the state of the last 4ti2 API call.
_4ti2_status _4ti2_errno = _4ti2_OK;

_4ti2_status
_4ti2_get_errno()
{
    return _4ti2_errno;
}

_4ti2_state*
_4ti2_qsolve_create_state(_4ti2_precision prec)
{
    _4ti2_errno = _4ti2_OK;
    return new _4ti2_::QSolveAPI();
}

_4ti2_state*
_4ti2_rays_create_state(_4ti2_precision prec)
{
    _4ti2_errno = _4ti2_OK;
    return new _4ti2_::QSolveAPI();
}

_4ti2_state*
_4ti2_circuits_create_state(_4ti2_precision prec)
{
    _4ti2_errno = _4ti2_OK;
    return new _4ti2_::QSolveAPI();
}

_4ti2_status
_4ti2_state_set_options(_4ti2_state* state, int argc, char** argv)
{
    state->set_options(argc, argv);
    return (_4ti2_errno = _4ti2_OK);
}

void
_4ti2_state_delete(_4ti2_state* state)
{
    _4ti2_errno = _4ti2_OK;
    delete state;
}

_4ti2_status
_4ti2_state_compute(_4ti2_state* state)
{
    state->compute();
    return (_4ti2_errno = _4ti2_OK);
}

_4ti2_matrix*
_4ti2_state_create_matrix(_4ti2_state* state, int num_rows, int num_cols, const char* name)
{
    _4ti2_errno = _4ti2_OK;
    return state->create_matrix(num_rows, num_cols, name);
}

_4ti2_matrix*
_4ti2_state_read_matrix(_4ti2_state* state, const char* filename, const char* name)
{
    _4ti2_errno = _4ti2_OK;
    return state->create_matrix(filename, name);
}

_4ti2_matrix*
_4ti2_state_get_matrix(_4ti2_state* state, const char* name)
{
    _4ti2_errno = _4ti2_OK;
    return state->get_matrix(name);
}

int
_4ti2_matrix_get_num_rows(_4ti2_matrix*  matrix)
{
    _4ti2_errno = _4ti2_OK;
    return matrix->get_num_rows();
}

int
_4ti2_matrix_get_num_cols(_4ti2_matrix*  matrix)
{
    _4ti2_errno = _4ti2_OK;
    return matrix->get_num_cols();
}

void
_4ti2_matrix_write_to_stdout(_4ti2_matrix*  matrix)
{
    matrix->write(std::cout);
    _4ti2_errno = _4ti2_OK;
}

void
_4ti2_matrix_write_to_stderr(_4ti2_matrix*  matrix)
{
    matrix->write(std::cerr);
    _4ti2_errno = _4ti2_OK;
}

void
_4ti2_matrix_write_to_file(_4ti2_matrix*  matrix, const char* filename)
{
    matrix->write(filename);
    _4ti2_errno = _4ti2_OK;
}

_4ti2_status
_4ti2_matrix_set_entry_int32_t(_4ti2_matrix*  matrix, int r, int c, int32_t value)
{
    matrix->set_entry_int32_t(r, c, value);
    return (_4ti2_errno = _4ti2_OK);
}

_4ti2_status
_4ti2_matrix_get_entry_int32_t(_4ti2_matrix*  matrix, int r, int c, int32_t* value)
{
    matrix->get_entry_int32_t(r, c, *value);
    return (_4ti2_errno = _4ti2_OK);
}

_4ti2_status
_4ti2_matrix_set_entry_int64_t(_4ti2_matrix*  matrix, int r, int c, int64_t value)
{
    matrix->set_entry_int64_t(r, c, value);
    return (_4ti2_errno = _4ti2_OK);
}

_4ti2_status
_4ti2_matrix_get_entry_int64_t(_4ti2_matrix*  matrix, int r, int c, int64_t* value)
{
    matrix->get_entry_int64_t(r, c, *value);
    return (_4ti2_errno = _4ti2_OK);
}

#ifdef _4ti2_GMP_
_4ti2_status
_4ti2_matrix_set_entry_mpz_ptr(_4ti2_matrix*  matrix, int r, int c, mpz_ptr value)
{
    mpz_class z(value);
    matrix->set_entry_mpz_class(r, c, z);
    return (_4ti2_errno = _4ti2_OK);
}

_4ti2_status
_4ti2_matrix_get_entry_mpz_ptr(_4ti2_matrix*  matrix, int r, int c, mpz_ptr value)
{
    mpz_class z;
    matrix->get_entry_mpz_class(r, c, z);
    mpz_set(value, z.get_mpz_t());
    return (_4ti2_errno = _4ti2_OK);
}
#endif

}

