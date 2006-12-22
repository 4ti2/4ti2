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

#include "FilterReduction.h"
#include "BinomialStream.h"
#include "Filter.h"
#include "Globals.h"

#include <vector>

#define END Binomial::get_num_svars()

using namespace _4ti2_;

typedef std::vector<const Binomial*> BinomialList;

class FilterNode {
public:
    FilterNode() { binomials = 0; filter = 0; }
    virtual ~FilterNode()
    {
        delete binomials; delete filter;
        for (int i = 0; i < (int) nodes.size(); ++i) { delete nodes[i].second; }
    }
    std::vector<std::pair<int,FilterNode*> > nodes;
    BinomialList* binomials;
    Filter* filter;
};

FilterReduction::FilterReduction()
{
    root = new FilterNode();
}

FilterReduction::~FilterReduction()
{
    delete root;
}

void
FilterReduction::add(const Binomial& b)
{
    FilterNode* current = root;
    Index end = END;
    for (Index i = 0; i < end; ++i)
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
                FilterNode* next = new FilterNode;
                current->nodes.push_back(std::pair<int,FilterNode*>(i,next));
                current = next;
            }
        }
    }
    if (!current->binomials)
    {
        current->binomials = new BinomialList;
        current->filter = new Filter;
        b.get_filter(*current->filter);
    }
    current->binomials->push_back(&b);
}

// Assumes point exists.
void
FilterReduction::remove(const Binomial& b)
{
    FilterNode* current = root;
    Index end = END;
    for (Index i = 0; i < end; ++i)
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
        if ((*i) == &b)
        {
            current->binomials->erase(i);
            return;
        }
    }
    assert(false); // If we got to here, then we did not find the point.
}

void
FilterReduction::clear()
{
    delete root;
    root = new FilterNode();
}

const Binomial*
FilterReduction::reducable(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable(b, b1, root);
}

const Binomial*
FilterReduction::reducable_negative(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable_negative(b, b1, root);
}

const Binomial*
FilterReduction::reducable(
                    const Binomial& b,
                    const Binomial* b1,
                    const FilterNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        if (b[node->nodes[i].first] > 0)
        {
            const Binomial* bi = reducable(b, b1, node->nodes[i].second);
            if (bi != 0) { return bi; }
        }
    }

    if (node->binomials)
    {
        Filter& f = *node->filter;
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            const Binomial& bi = *(*i);
            if (Binomial::reduces(bi, f, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

const Binomial*
FilterReduction::reducable_negative(
                    const Binomial& b,
                    const Binomial* b1,
                    const FilterNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        if (b[node->nodes[i].first] < 0)
        {
            const Binomial* bi = reducable_negative(b, b1, node->nodes[i].second);
            if (bi != 0) { return bi; }
        }
    }

    if (node->binomials)
    {
        Filter& f = *node->filter;
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            const Binomial& bi = *(*i);
            if (Binomial::reduces_negative(bi, f, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

void
FilterReduction::print() const
{
    print(root);
}

void
FilterReduction::print(const FilterNode* node) const
{
    assert(node != 0);
    if (node->binomials)
    {
        *out << "Num binomials = " << node->binomials->size() << std::endl;
        for (int i = 0; i < (int) node->filter->size(); ++i)
        {
            *out << (*node->filter)[i] << " ";
        }
        *out << "\n";
        for (BinomialList::iterator i = node->binomials->begin();
              i != node->binomials->end(); ++i)
        {
            *out << *(*i) << "\n";
        }
    }

    for (int i = 0; i < (int) node->nodes.size(); ++i)
    {
        print(node->nodes[i].second);
    }
}

