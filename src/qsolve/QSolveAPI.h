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

#ifndef _4ti2_qsolve__QSolveAPI_
#define _4ti2_qsolve__QSolveAPI_

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"
#include "qsolve/QSolveVariant.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/VectorArray.h"

namespace _4ti2_ {

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
    std::string basename;
    _4ti2_precision prec;
    QSolveVariant algorithm;
    QSolveConsOrder order;
    _4ti2_constraint sign_default;
    _4ti2_constraint rel_default;
    enum _4ti2_input { DEFAULT, INE, EXT, IEQ, POI };
    _4ti2_input input;

    void initialise_data();
    template <class T>
    void initialise_dataT();
    template <class T>
    void computeT();

    virtual void pre_compute();
    virtual void post_compute();

    virtual void write_usage();
    virtual void write_options();
    virtual void write_input_files();
    virtual void write_output_files();

    void parse_porta_ieq(std::istream&in);

    void print_banner();
    void unrecognised_option_argument(const char* option);

    _4ti2_matrix* mat;
    _4ti2_matrix* sign;
    _4ti2_matrix* rel;
    _4ti2_matrix* ray;
    _4ti2_matrix* cir;
    _4ti2_matrix* qhom;
    _4ti2_matrix* qfree;
};

} // namespace _4ti2_

#endif
