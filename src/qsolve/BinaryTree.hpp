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

////////////////
// BinaryTree //
////////////////

template <class IndexSet>
const int BinaryTree<IndexSet>::INDENT = 2;

template <class IndexSet>
BinaryTree<IndexSet>::BinaryTree()
        : supp_size(0)
{
    root = new TreeBucket();
}

template <class IndexSet>
BinaryTree<IndexSet>::~BinaryTree()
{
    clear();
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::clear()
{
    delete root;
    root = new TreeBucket();
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::insert(const std::vector<IndexSet>& supports)
{
    insert(supports, 0, supports.size());
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::insert(const std::vector<IndexSet>& supports, Index start, Index end)
{
    for (int i = start; i < end; ++i) { insert(supports[i], i); }
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::insert(const IndexSet& support, Index index)
{
    supp_size = support.get_size();
    root = root->insert(support, index);
}

template <class IndexSet>
inline
bool
BinaryTree<IndexSet>::dominated(const IndexSet& support, Index index1, Index index2) const
{
    return root->dominated(support, index1, index2);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::dump() const
{
    *out << "Binary Tree Dump:\n";
    root->dump(0);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    root->find_singleton_diff(inds, s);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    root->find(inds, zeros, s, count);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::find(std::vector<std::pair<Index,Index> >& inds, const BinaryTree& tree, Size count) const
{
    IndexSet s(supp_size,0);
    root->find(inds, s, tree.root, count);
}

#if 0
////////////////////
// BinaryTreeLeaf //
////////////////////

template <class IndexSet>
typename BinaryTree<IndexSet>::TreeNode*
BinaryTree<IndexSet>::BinaryTreeLeaf::insert(const IndexSet& s, Index index)
{
    IndexSet tmp(s);
    tmp.set_symm_difference(is);
    Index new_index = 0;
    while (!tmp[new_index]) { ++new_index; }

    BinaryTreeLeaf* new_leaf = new BinaryTreeLeaf(index, s);
    TreeBranch* new_branch;
    if (s[new_index]) {
        new_branch = new TreeBranch(new_index, is, this, s, new_leaf);
    }
    else {
        new_branch = new TreeBranch(new_index, s, new_leaf, is, this);
    }
    return new_branch;
}

template <class IndexSet>
bool
BinaryTree<IndexSet>::BinaryTreeLeaf::dominated(const IndexSet& b, Index index1, Index index2) const
{
    if (i != index1 && i != index2) { return true; }
    return false;
}

template <class IndexSet>
void
BinaryTree<IndexSet>::BinaryTreeLeaf::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "[" << i << "] " << is << "\n";
}
#endif

//////////////////////
// TreeBranch //
//////////////////////

template <class IndexSet>
BinaryTree<IndexSet>::TreeBranch::TreeBranch(Index index, const IndexSet& supp0, TreeNode* _zero, 
                const IndexSet& supp1, TreeNode* _one)
            : i(index), is0(supp0), is1(supp1), zero(_zero), one(_one)
{
}

template <class IndexSet>
BinaryTree<IndexSet>::TreeBranch::~TreeBranch()
{
    delete one;
    delete zero;
}

template <class IndexSet>
typename BinaryTree<IndexSet>::TreeNode*
BinaryTree<IndexSet>::TreeBranch::insert(const IndexSet& s, Index index)
{
    if (s[i]) { is1.set_intersection(s); one = one->insert(s, index); }
    else { is0.set_intersection(s); zero = zero->insert(s, index); }
    return this;
}

template <class IndexSet>
bool
BinaryTree<IndexSet>::TreeBranch::dominated(const IndexSet& b, Index index1, Index index2) const
{
    if (is0.set_subset(b) && zero->dominated(b, index1, index2)) { return true; }
    if (is1.set_subset(b) && one->dominated(b, index1, index2)) { return true; }
    return false;
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBranch::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    if (is0.singleton_diff(s)) { zero->find_singleton_diff(inds, s); }
    if (is1.singleton_diff(s)) { one->find_singleton_diff(inds, s); }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBranch::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    if (is0.set_disjoint(zeros) && is0.count_union(s) <= count) { zero->find(inds, zeros, s, count); }
    if (is1.set_disjoint(zeros) && is1.count_union(s) <= count) { one->find(inds, zeros, s, count); }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBranch::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, const TreeNode* tree, Size count) const
{
    if (is0.count_union(supp) <= count) { tree->find(inds, is0, zero, count); }
    if (is1.count_union(supp) <= count) { tree->find(inds, is1, one, count); }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBranch::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, Index index, Size count) const
{
    if (is0.count_union(supp) <= count) { zero->find(inds, supp, index, count); }
    if (is1.count_union(supp) <= count) { one->find(inds, supp, index, count); }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBranch::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << i << " = 0:\n";
    zero->dump(level+1);
    *out << indent << i << " = 1:\n";
    one->dump(level+1);
}

//////////////////////
// TreeBucket //
//////////////////////

template <class IndexSet>
const Size BinaryTree<IndexSet>::TreeBucket::MAX_BUCKET_SIZE = 50;

template <class IndexSet>
BinaryTree<IndexSet>::TreeBucket::TreeBucket()
{
}

template <class IndexSet>
BinaryTree<IndexSet>::TreeBucket::~TreeBucket()
{
}

template <class IndexSet>
typename BinaryTree<IndexSet>::TreeNode*
BinaryTree<IndexSet>::TreeBucket::insert(const IndexSet& s, Index index)
{
    supps.push_back(s);
    indices.push_back(index);

    if ((Size) supps.size() >= MAX_BUCKET_SIZE) {
#if 0
        IndexSet tmp(supps[0]);
        tmp.set_symm_difference(supps[1]);
        Index new_index = 0;
        while (!tmp[new_index]) { ++new_index; }
#endif

        int new_count = supps.size()+1;
        Index new_index = -1;
        for (Index i = 0; i < s.get_size(); ++i) {
            int count = -supps.size()/2;
            for (Index j = 0; j < (Index) supps.size(); ++j) { if (supps[j][i]) { ++count; } }
            if (count < 0) { count *= -1; }
            if (count < new_count) { new_index = i; new_count = count; }
        }

        IndexSet is0(s.get_size(),1);
        IndexSet is1(s.get_size(),1);
        TreeBucket* one_bucket = new TreeBucket();
        for (Index i = supps.size()-1; i >= 0; --i) {
            if (supps[i][new_index]) {
                one_bucket->insert(supps[i], indices[i]);
                is1.set_intersection(supps[i]);
                supps.erase(supps.begin()+i);
                indices.erase(indices.begin()+i);
            }
            else {
                is0.set_intersection(supps[i]);
            }
        }
        TreeBranch* new_branch = new TreeBranch(new_index, is0, this, is1, one_bucket);
        return new_branch;
    }
    return this;
}

template <class IndexSet>
bool
BinaryTree<IndexSet>::TreeBucket::dominated(const IndexSet& s, Index index1, Index index2) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_subset(s) && indices[i] != index1 && indices[i] != index2) { return true; }
    }
    return false;
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBucket::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].singleton_diff(s)) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBucket::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_disjoint(zeros) && supps[i].count_union(s) <= count) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBucket::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, const TreeNode* tree, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].count_union(s) <= count) { tree->find(inds, supps[i], indices[i], count); }
    }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBucket::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, Index index, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].count_union(s) <= count) { inds.push_back(std::pair<Index,Index>(indices[i],index)); }
    }
}

template <class IndexSet>
void
BinaryTree<IndexSet>::TreeBucket::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        *out << indent << "[" << indices[i] << "] " << supps[i] << "\n";
    }
}
