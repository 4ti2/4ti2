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

#include "qsolve/BinaryTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

template <class IndexSet>
const int BinaryTree<IndexSet>::INDENT = 2;

template <class IndexSet>
BinaryTree<IndexSet>::BinaryTreeBranch::BinaryTreeBranch(Index index, BinaryTreeNode* _zero, BinaryTreeNode* _one)
            : i(index), zero(_zero), one(_one)
{
}

template <class IndexSet>
BinaryTree<IndexSet>::BinaryTreeBranch::~BinaryTreeBranch()
{
    delete one;
    delete zero;
}

template <class IndexSet>
typename BinaryTree<IndexSet>::BinaryTreeNode*
BinaryTree<IndexSet>::BinaryTreeBranch::insert(const IndexSet& support, int index)
{
    if (support[i]) { one = one->insert(support, index); }
    else { zero = zero->insert(support, index); }
    return this;
}

template <class IndexSet>
typename BinaryTree<IndexSet>::BinaryTreeNode*
BinaryTree<IndexSet>::BinaryTreeLeaf::insert(const IndexSet& support, int index)
{
    IndexSet tmp(support);
    tmp.set_symm_difference(is);
    Index new_index = 0;
    while (!tmp[new_index]) { ++new_index; }

    BinaryTreeLeaf* new_leaf = new BinaryTreeLeaf(index, support);
    BinaryTreeBranch* new_branch;
    if (support[new_index]) {
        new_branch = new BinaryTreeBranch(new_index, this, new_leaf);
    }
    else {
        new_branch = new BinaryTreeBranch(new_index, new_leaf, this);
    }
    return new_branch;
}

template <class IndexSet>
bool
BinaryTree<IndexSet>::BinaryTreeBranch::dominated(const IndexSet& b, int index1, int index2)
{
    if (zero->dominated(b, index1, index2)) { return true; }
    if (b[i] && one->dominated(b, index1, index2)) { return true; } 
    return false;
}

template <class IndexSet>
void
BinaryTree<IndexSet>::BinaryTreeBranch::find_diff(
            std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff)
{
    zero->find_diff(indices, supp, diff, orig_diff);
    if (supp[i]) {
        if (diff > 0) { one->find_diff(indices, supp, diff-1, orig_diff); }
    }
    else {
        one->find_diff(indices, supp, diff, orig_diff);
    }
}
