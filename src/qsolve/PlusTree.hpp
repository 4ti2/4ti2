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

#include "qsolve/PlusTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

template <class IndexSet>
const int PlusTree<IndexSet>::INDENT = 2;

template <class IndexSet>
PlusTree<IndexSet>::TreeBranch::TreeBranch()
{
}

template <class IndexSet>
PlusTree<IndexSet>::TreeBranch::~TreeBranch()
{
    clear();
}

template <class IndexSet>
inline void 
PlusTree<IndexSet>::TreeBranch::clear()
{
    for (unsigned int i = 0; i < nodes.size(); ++i)  { delete nodes[i].second; }
    nodes.clear();
}

template <class IndexSet>
typename PlusTree<IndexSet>::TreeNode*
PlusTree<IndexSet>::TreeBranch::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);

    // There should always be another.
    int next_one = start;
    while (!support[next_one]) { ++next_one; }

    size_t i = 0;
    while (i < nodes.size() && next_one > nodes[i].first) { ++i; }
    if (i < nodes.size() && next_one == nodes[i].first) {
        nodes[i].second = nodes[i].second->insert(support, next_one+1, index);
    }
    else {
        TreeNode* new_node = new TreeLeaf(index, support);
        nodes.insert(nodes.begin()+i, std::pair<int,TreeNode*>(next_one,new_node));
    }
    return this;

#if 0
    unsigned i = 0;
    while (i < nodes.size() && next_one != nodes[i].first) { ++i; }
    if (i < nodes.size()) {
        nodes[i].second = nodes[i].second->insert(support, next_one+1, index);
    }
    else {
        TreeLeaf* new_node = new TreeLeaf(index, support);
        nodes.push_back(std::pair<int,TreeNode*>(next_one,new_node));
    }
    return this;
#endif
}

template <class IndexSet>
typename PlusTree<IndexSet>::TreeNode*
PlusTree<IndexSet>::TreeLeaf::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);

    // There should always be another.
    int next_one = start;
    while (!is[next_one]) { ++next_one; }

    TreeBranch* new_branch = new TreeBranch;
    new_branch->add(next_one, this);
    new_branch->insert(support, start, index);

    return new_branch;
}

template <class IndexSet>
bool
PlusTree<IndexSet>::TreeBranch::dominated(const IndexSet& b, int index1, int index2) const
{
    for (unsigned i = 0; i < nodes.size(); ++i) {
        if (b[nodes[i].first] && nodes[i].second->dominated(b, index1, index2)) {
            return true;
        }
    }
    return false;
}

template <class IndexSet>
void
PlusTree<IndexSet>::TreeBranch::find_diff(
            std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff) const
{
    for (unsigned i = 0; i < nodes.size(); ++i) {
        if (supp[nodes[i].first]) {
            if (diff > 0) { nodes[i].second->find_diff(indices, supp, diff-1, orig_diff); }
        }
        else {
            nodes[i].second->find_diff(indices, supp, diff, orig_diff);
        }
    }
}
