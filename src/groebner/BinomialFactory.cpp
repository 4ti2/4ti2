/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
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

#include "BinomialFactory.h"
#include "WeightAlgorithm.h"
#include <iostream>
#include "BitSetStream.h"
#include "VectorArrayStream.h"
#include "VectorStream.h"
#include "Globals.h"

#include "Debug.h"

using namespace _4ti2_;

BinomialFactory::BinomialFactory(
                Feasible& _feasible,
                const VectorArray& cost)
{
    permutation = 0;
    costs = 0;
    orig_bnd = 0;

    VectorArray tmp_cost(cost);
    check_cost(_feasible, tmp_cost);

    initialise( _feasible.get_dimension(),
                _feasible.get_basis(),
                tmp_cost,
                _feasible.get_urs(),
                _feasible.get_bnd(),
                _feasible.get_unbnd(),
                _feasible.get_grading(),
                _feasible.get_weights(),
                _feasible.get_max_weights(),
                _feasible.get_rhs());
}

BinomialFactory::BinomialFactory(
                Feasible& _feasible,
                const VectorArray& cost,
                const BitSet& sat)
{
    permutation = 0;
    costs = 0;
    orig_bnd = 0;

    VectorArray tmp_cost(cost);
    check_cost(_feasible, tmp_cost);

    initialise( _feasible.get_dimension(),
                _feasible.get_basis(),
                tmp_cost,
                _feasible.get_urs(),
                sat,
                _feasible.get_unbnd(),
                _feasible.get_grading(),
                _feasible.get_weights(),
                _feasible.get_max_weights(),
                _feasible.get_rhs());
}

// Checks whether the cost function is bounded.
// It adds a vector to cost if needed to ensure that it is bounded.
void
BinomialFactory::check_cost(Feasible feasible, VectorArray& cost)
{
    BitSet cost_unbnd(feasible.get_dimension());
    bool is_bounded = feasible.bounded(cost, cost_unbnd);
    if (!is_bounded)
    {
        std::cerr << "Cost function is not bounded.\n";
        exit(1);
    }
    else if (!cost_unbnd.empty())
    {
        // If there are unbounded components, then we must insert another cost
        // vector to ensure that the miniumum is bounded.
        Vector tmp(cost.get_size(),0);
        for (int i = 0; i < cost.get_size(); ++i)
        {   if (cost_unbnd[i]) { tmp[i] = 1; } }
        cost.insert(tmp);
    }
}

BinomialFactory::~BinomialFactory()
{
    delete permutation;
    delete costs;
    delete orig_bnd;
}

void
BinomialFactory::initialise(
                int dim,
                const VectorArray& _lattice,
                const VectorArray& cost,
                const BitSet& urs,
                const BitSet& bnd,
                const BitSet& unbnd,
                const Vector& _grading,
                const VectorArray* _weights,
                const Vector* _max_weights,
                const Vector* _rhs)
{
    assert(urs.get_size() == dim);
    assert(bnd.get_size() == dim);
    assert(_grading.get_size() == dim);
    assert(cost.get_size() == dim);

    delete orig_bnd;
    orig_bnd = new BitSet(bnd);

    delete costs;
    costs = new VectorArray(cost);

    Binomial::bnd_end = bnd.count();
    Binomial::rs_end = dim - urs.count();
    Binomial::size = dim + costs->get_number();
    Binomial::urs_end = dim;
    Binomial::cost_start = dim;
    Binomial::cost_end = Binomial::size;

    delete permutation;
    initialise_permutation(bnd, urs);

    delete Binomial::grading;
    Binomial::grading = new Grading(_grading);
    Binomial::grading->permute(*permutation);
    DEBUG_4ti2(std::cout << "Grading:\n" << *Binomial::grading << "\n";)

    set_weights(_weights, _max_weights);
    set_truncated(_lattice, _rhs);
}

void
BinomialFactory::set_truncated(
                const VectorArray& _lattice,
                const Vector* _rhs)
{
    delete Binomial::rhs; Binomial::rhs = 0;
    delete Binomial::lattice; Binomial::lattice = 0;
    if (Globals::truncation == Globals::NONE) { return; }
    if (_rhs != 0 && orig_bnd->count() != 0)
    {
        if (Globals::truncation != Globals::WEIGHT)
        {
            assert(_rhs->get_size() == _lattice.get_size());
            assert(_rhs->get_size() == Binomial::urs_end);
            Binomial::rhs = new Vector(orig_bnd->count());
            Vector::project(*_rhs, *orig_bnd, *Binomial::rhs);
            Binomial::lattice = new VectorArray(_lattice.get_number(),
                            orig_bnd->count());
            VectorArray::project(_lattice, *orig_bnd, *Binomial::lattice);
            DEBUG_4ti2(*out<<"Truncation Lattice:\n"<<*Binomial::lattice<<"\n";)
            DEBUG_4ti2(*out<<"Truncation RHS:\n"<<*Binomial::rhs<<"\n";)
        }

        BitSet unbnd(*orig_bnd);
        unbnd.set_complement();
        Vector weight(_lattice.get_size(),0);
        Vector zero(_lattice.get_size(),0);
        if (Globals::norm == 2)
        {
            lp_weight_l2(_lattice, unbnd, *_rhs, weight);
        }
        else
        {
            lp_weight_l1(_lattice, unbnd, *_rhs, weight);
        }
        IntegerType max = Vector::dot(*_rhs, weight);
        DEBUG_4ti2(*out<<"Weight:\n"<<weight<<"\n";)
        if (weight != zero) { add_weight(weight, max); }
    }
}

void
BinomialFactory::set_weights(
                const VectorArray* _weights,
                const Vector* _max_weights)
{
    delete Binomial::weights; Binomial::weights = 0;
    delete Binomial::max_weights; Binomial::max_weights = 0;
    if (_weights != 0 && _max_weights != 0)
    {
        assert(_weights->get_number() == _max_weights->get_size());
        assert(_weights->get_size() == Binomial::urs_end);
        Binomial::weights = new VectorArray(*_weights);
        Binomial::max_weights = new Weight(*_max_weights);
        BitSet unbnd(*orig_bnd);
        unbnd.set_complement();
        WeightAlgorithm::strip_weights(Binomial::weights, Binomial::max_weights, unbnd);
        Binomial::weights->permute(*permutation);
        DEBUG_4ti2(*out<<"Weights:\n"<<*Binomial::weights<<"\n";)
        DEBUG_4ti2(*out<<"RHS:\n"<<*Binomial::max_weights<<"\n";)
    }
}

void
BinomialFactory::add_weight(
                const Vector& _weight,
                IntegerType _max_weight)
{
    Vector tmp(_weight);
    tmp.permute(*permutation);
    DEBUG_4ti2(*out<<"Weight:\n"<<tmp<<"\n";)
    if (Binomial::weights == 0 || Binomial::max_weights == 0)
    {
        Binomial::weights = new VectorArray(0, _weight.get_size());
        Binomial::weights->insert(tmp);
        Binomial::max_weights = new Vector(1, _max_weight);
    }
    else
    {
        Binomial::weights->insert(tmp);
        Vector tmp_max(1, _max_weight);
        Vector* new_max_weights = new Vector(Binomial::max_weights->get_size()+1);
        Vector::concat(*Binomial::max_weights, tmp_max, *new_max_weights);
        delete Binomial::max_weights;
        Binomial::max_weights = new_max_weights;
    }
}

void
BinomialFactory::initialise_permutation(
                const BitSet& bnd_mask,
                const BitSet& urs_mask)
{
    // Check that no variable is both urs and bounded.
    assert(BitSet::set_disjoint(bnd_mask,urs_mask));

    int num_bnd = bnd_mask.count();
    int num_urs = urs_mask.count();

    permutation = new Permutation(bnd_mask.get_size());

    int bnd_count = 0;
    int urs_count = 0;
    int nothing_count = 0;
    // TODO: Change name of locals.
    int bnd_index = 0;
    int urs_index = bnd_mask.get_size() - num_urs;
    int nothing_index = num_bnd;
    for (Index i = 0; i < bnd_mask.get_size(); ++i)
    {
        if (urs_mask[i] == 1)
        {
            (*permutation)[urs_index + urs_count] = i;
            ++urs_count;
        }
        else if (bnd_mask[i] == 1)
        {
            (*permutation)[bnd_index + bnd_count] = i;
            ++bnd_count;
        }
        else
        {
            (*permutation)[nothing_index + nothing_count] = i;
            ++nothing_count;
        }
    }

    // TODO: Check to see if a permutation is necessary.
    // i.e. Is the permutation trivial?
}

void
BinomialFactory::convert(const Vector& v, Binomial& b) const
{
    assert(b.size == v.get_size() + costs->get_number());
    assert(v.get_size() == costs->get_size());

    for (Index i = 0; i < v.get_size(); ++i)
    {
        b[i] = v[(*permutation)[i]];
    }
    const VectorArray& cost_matrix = *costs;
    for (Index i = 0; i < cost_matrix.get_number(); ++i)
    {
        b[i+Binomial::cost_start] = Vector::dot(v,cost_matrix[i]);
    }
}

void
BinomialFactory::convert(const Binomial& b, Vector& v) const
{
    assert(v.get_size() == b.urs_end);
    for (Index i = 0; i < v.get_size(); ++i)
    {
        v[(*permutation)[i]] = b[i];
    }
}

void
BinomialFactory::convert(const VectorArray& vs, BinomialCollection& bc, bool orientate) const
{
    Binomial b;
    for (int i = 0; i < vs.get_number(); ++i)
    {
        convert(vs[i], b);
        if (!Binomial::overweight(b) && !Binomial::truncated(b))
        {
            if (orientate) { if (b.orientate()) { bc.add(b); } }
            else { bc.add(b); }
            assert(!b.is_non_positive());
        }
    }
}

void
BinomialFactory::convert(const BinomialArray& bs, VectorArray& vs) const
{
    assert(vs.get_size() == Binomial::urs_end);
    vs.renumber(bs.get_number());
    for (int i = 0; i < bs.get_number(); ++i) { convert(bs[i], vs[i]); }
}

void
BinomialFactory::convert(const BinomialSet& bs, VectorArray& vs) const
{
    assert(vs.get_size() == Binomial::urs_end);
    vs.renumber(bs.get_number());
    for (int i = 0; i < bs.get_number(); ++i) { convert(bs[i], vs[i]); }
}
