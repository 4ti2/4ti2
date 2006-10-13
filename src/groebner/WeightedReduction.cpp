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

#include "WeightedReduction.h"
#include "BinomialStream.h"
#include "Globals.h"

#include <vector>
#include <map>

#define END Binomial::get_num_svars()

using namespace _4ti2_;

typedef std::multimap<IntegerType, const Binomial*> BinomialList;

class WeightedNode {
public:
    WeightedNode() { binomials = 0; }
    virtual ~WeightedNode() { delete binomials; }
    std::vector<std::pair<int,WeightedNode*> > nodes;
    BinomialList* binomials;
};

WeightedReduction::WeightedReduction()
{
    root = new WeightedNode();
}

WeightedReduction::~WeightedReduction()
{
    delete root;
}

void
WeightedReduction::add(const Binomial& b)
{
    WeightedNode* current = root;
    Index end = END;
    for (Index i = 0; i < end-1; ++i)
    {
        if (b[i] > 0)
        {
            int j = 0;
            while (j < (int) current->nodes.size() && current->nodes[j].first != i) { ++j; }
            if (j < (int) current->nodes.size())
            {
                current = current->nodes[j].second;
            }
            else
            {
                WeightedNode* next = new WeightedNode;
                current->nodes.push_back(std::pair<int,WeightedNode*>(i,next));
                current = next;
            }
        }
    }
    if (!current->binomials) { current->binomials = new BinomialList; }
    current->binomials->insert(BinomialList::value_type(b.l1_norm(), &b));
}

// Assumes point exists.
void
WeightedReduction::remove(const Binomial& b)
{
    WeightedNode* current = root;
    Index end = END;
    for (Index i = 0; i < end-1; ++i)
    {
        if (b[i] > 0)
        {
            int j = 0;
            while (j < (int) current->nodes.size() && current->nodes[j].first != i) { ++j; }
            if (j < (int) current->nodes.size())
            {
                current = current->nodes[j].second;
            }
            else
            {
                assert(false); // If we got to here, then we did not find the point.
            }
        }
    }
    assert(current->binomials); // If this assert is false, the point does not exist.
    for (BinomialList::iterator i = current->binomials->begin();
                    i != current->binomials->end(); ++i)
    {
        if ((*i).second == &b)
        {
            current->binomials->erase(i);
            return;
        }
    }
    assert(false); // If we got to here, then we did not find the point.
}

void
WeightedReduction::clear()
{
    delete root;
    root = new WeightedNode();
}

const Binomial*
WeightedReduction::reducable(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable(b, b.l1_norm(), b1, root);
}

const Binomial*
WeightedReduction::reducable_negative(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable_negative(b, b.l1_norm_negative(), b1, root);
}

const Binomial*
WeightedReduction::reducable(
                    const Binomial& b,
                    const IntegerType& norm,
                    const Binomial* b1,
                    const WeightedNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        if (b[node->nodes[i].first] > 0)
        {
            const Binomial* bi = reducable(b, norm, b1, node->nodes[i].second);
            if (bi != 0) { return bi; }
        }
    }

    if (node->binomials)
    {
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            if (norm < (*i).first) { break; }
            const Binomial& bi = *(*i).second;
            if (Binomial::reduces(bi, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

const Binomial*
WeightedReduction::reducable_negative(
                    const Binomial& b,
                    const IntegerType& norm,
                    const Binomial* b1,
                    const WeightedNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        if (b[node->nodes[i].first] < 0)
        {
            const Binomial* bi = reducable_negative(b, norm, b1, node->nodes[i].second);
            if (bi != 0) { return bi; }
        }
    }

    if (node->binomials)
    {
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            if (norm < (*i).first) { break; }
            const Binomial& bi = *(*i).second;
            if (Binomial::reduces_negative(bi, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

void
WeightedReduction::print() const
{
    print(root);
}

void
WeightedReduction::print(const WeightedNode* node) const
{
    assert(node != 0);
    if (node->binomials)
    {
        *out << "Num binomials = " << node->binomials->size() << std::endl;
        for (BinomialList::iterator i = node->binomials->begin();
              i != node->binomials->end(); ++i)
        {
            *out << *(*i).second << "\n";
        }
    }

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        print(node->nodes[i].second);
    }
}

