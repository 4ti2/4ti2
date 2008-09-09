#ifndef _4ti2API_
#define _4ti2API_

#include <iostream>
#ifdef _4ti2_GMP_
#include <gmp.h>
#include <gmpxx.h>
#endif

struct _4ti2_matrix {
public:
    _4ti2_matrix() {}
    virtual ~_4ti2_matrix() {}

    virtual int get_num_rows() const = 0;
    virtual int get_num_cols() const = 0;

    virtual void write(const char* filename) const = 0;
    virtual void write(std::ostream& out) const = 0; 
    virtual void read(std::istream& in) = 0; 

    virtual void set_entry_int32_t(int r, int c, const int32_t& value) = 0; 
    virtual void get_entry_int32_t(int r, int c, int32_t& value) const = 0;
    virtual void set_entry_int64_t(int r, int c, const int64_t& value) = 0;
    virtual void get_entry_int64_t(int r, int c, int64_t& value) const = 0;

#ifdef _4ti2_GMP_
    virtual void set_entry_mpz_class(int r, int c, const mpz_class& value) = 0;
    virtual void get_entry_mpz_class(int r, int c, mpz_class& value) const = 0;
#endif
};

struct _4ti2_state {
public:
    _4ti2_state() {}
    virtual ~_4ti2_state() {}

    virtual void compute() = 0;

    virtual void set_options(int argc, char** argv) = 0; 

    virtual void read(const char* project) = 0;
    virtual void write(const char* project) = 0;

    virtual _4ti2_matrix* create_matrix(int num_rows, int num_cols, const char* name) = 0;
    virtual _4ti2_matrix* create_matrix(const char* filename, const char* name) = 0;
    virtual _4ti2_matrix* create_matrix(std::istream& in, const char* name) = 0;

    virtual _4ti2_matrix* get_matrix(const char* name) = 0;
};

#endif
