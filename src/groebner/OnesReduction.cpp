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

#include "OnesReduction.h"
#include "BinomialStream.h"
#include "Globals.h"

#include <vector>

#define END Binomial::get_num_svars()

using namespace _4ti2_;

typedef std::vector<const Binomial*> BinomialList;

class OnesNode {
public:
    OnesNode() { binomials = 0; }
    virtual ~OnesNode() { delete binomials; }
    std::vector<std::pair<int,OnesNode*> > nodes;
    BinomialList* binomials;
};

OnesReduction::OnesReduction()
{
    root = new OnesNode();
}

OnesReduction::~OnesReduction()
{
    delete root;
}

void
OnesReduction::add(const Binomial& b)
{
    OnesNode* current = root;
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
                OnesNode* next = new OnesNode;
                current->nodes.push_back(std::pair<int,OnesNode*>(i,next));
                current = next;
            }
        }
    }
    if (!current->binomials) { current->binomials = new BinomialList; }
    current->binomials->push_back(&b);
}

// Assumes point exists.
void
OnesReduction::remove(const Binomial& b)
{
    OnesNode* current = root;
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
                // If we got to here, then we did not find the point.
                assert(false);
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
OnesReduction::clear()
{
    delete root;
    root = new OnesNode();
}

const Binomial*
OnesReduction::reducable(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable(b, b1, root);
}

const Binomial*
OnesReduction::reducable_negative(
                    const Binomial& b,
                    const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable_negative(b, b1, root);
}

const Binomial*
OnesReduction::reducable(
                    const Binomial& b,
                    const Binomial* b1,
                    const OnesNode* node) const
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
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            const Binomial& bi = *(*i);
            if (Binomial::reduces(bi, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

const Binomial*
OnesReduction::reducable_negative(
                    const Binomial& b,
                    const Binomial* b1,
                    const OnesNode* node) const
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
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i)
        {
            const Binomial& bi = *(*i);
            if (Binomial::reduces_negative(bi, b))
            {
                if (&bi != &b && &bi != b1) { return &bi; }
            }
        }
    }

    return 0;
}

void
OnesReduction::print() const
{
    print(root);
}

void
OnesReduction::print(const OnesNode* node) const
{
    assert(node != 0);
    if (node->binomials)
    {
        *out << "Num binomials = " << node->binomials->size() << std::endl;
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

