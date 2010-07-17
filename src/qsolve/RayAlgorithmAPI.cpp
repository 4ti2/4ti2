/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "qsolve/RayAlgorithmAPI.h"
#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/Cone.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

RayAlgorithmAPI::RayAlgorithmAPI()
{
}

void
RayAlgorithmAPI::compute()
{
    Timer t;
    *out << "Ray Matrix Algorithm.\n";
    DEBUG_4ti2(*out << "STATE:\n" << state_api << "\n";)

    // The number of variables.
    Size n = state_api->num_vars();
    // The number of constraints.
    Size m = state_api->num_cons();
    // The dimension.
    Size dim = state_api->num_gens();

    // The set of constraints to be processed.
    IndexSetD rem(state_api->num_vars()+state_api->num_cons(), 0);
    cone.constraint_set(_4ti2_LB, rem);
    for (Index i = 0; i < dim; ++i) { rem.unset(ineqs[i]); }

    // Construct the initial set of rays and supports.
    state_api->init();

    DEBUG_4ti2(*out << "Initial Dual:\n" << dual << "\n";)
    DEBUG_4ti2(*out << "Initial Supp:\n" << supp_api << "\n";)
    DEBUG_4ti2(*out << "Initial Rem:\n" << rem << "\n";)

    // The total set of relaxed constraints.
    IndexSetD rel(rem);
    cone.constraint_set(_4ti2_FR, rel);
    cone.constraint_set(_4ti2_DB, rel);

    Index next = -1;
    // Construct main algorithm object.
    RaySubAlgorithm* alg = construct_sub_algorithm(rel, cons_added, next);

    // Construct threaded algorithm objects.
    std::vector<RaySubAlgorithm*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg->clone()); }

    // While the cone is not empty and there are still rows to choose from.
    while (state_api->num_gens() > 0 && !rem.empty()) {
        DEBUG_4ti2(*out << "CONE:\n" << cone << "\n";)
        DEBUG_4ti2(*out << "DUAL:\n" << dual << "\n";)
        DEBUG_4ti2(*out << "SUPP:\n" << supps << "\n";)

        // Choose the next constraint and sort gens.
        Index pos_start, pos_end, neg_start, neg_end;
        next = state_api->next_constraint(pos_start, pos_end, neg_start, neg_end);

        // Check to see whether the next constraint is redundant. If so, ignore it.
        if (neg_start == neg_end) { rem.unset(next); continue; }

        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d", ENDL, rem.count(), next, rays.get_number());
        *out << buffer << std::flush;

        // Check to see if we need to perform any processing.
        if (pos_start != pos_end) {
            // Next, we assign index ranges for the inner and outer for loops so
            // that the outer loop is smaller than the inner loop.
            Index r1_start, r1_end, r2_start, r2_end;
            if (pos_end-pos_start <= neg_end-neg_start) {
                r1_start = pos_start; r1_end = pos_end;
                r2_start = neg_start; r2_end = neg_end;
            }
            else {
                r2_start = pos_start; r2_end = pos_end;
                r1_start = neg_start; r1_end = neg_end;
            }

            Index r2_index;
            // We sort the r2's into gens where r2_supp.count()==cons_added+1.
            state_api->sort_count(cons_added+1, r2_start, r2_end, r2_index);

            // Initialise master algorithm.
            alg.init(r1_start, r1_end, r2_start, r2_index, r2_end);
            // Run threads.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            // Run master algorithm.
            alg.compute();
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

            // Update the support vectors for the next_con.
            state_api->update(next, pos_start, pos_end);
            // Delete all the vectors with a negative entry in the column next.
            state_api->remove(neg_start, neg_end);

            // Add newly generated rays and supps to the dual cone and the supports.
            alg.transfer();
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }
    
        }

        rel.unset(next);
        rem.unset(next);
        ++cons_added;

        // Output statistics.
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                dual.num_gens(), t.get_elapsed_time());
        *out << buffer << std::endl;

        DEBUG_4ti2(*out << "CONE:\n" << cone << "\n";)
        DEBUG_4ti2(*out << "DUAL:\n" << dual << "\n";)
        DEBUG_4ti2(*out << "SUPP:\n" << supps << "\n";)
    }

    // Clean up algorithms objects.
    delete alg;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
}


