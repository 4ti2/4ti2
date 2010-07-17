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

#include "qsolve/MultiTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

////////////////
// MultiTree //
////////////////

template <class IndexSet>
const int MultiTree<IndexSet>::INDENT = 2;

template <class IndexSet>
MultiTree<IndexSet>::MultiTree()
        : supp_size(0)
{
    root = new TreeBucket();
}

template <class IndexSet>
MultiTree<IndexSet>::~MultiTree()
{
    clear();
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::clear()
{
    delete root;
    root = new TreeBucket();
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::insert(const std::vector<IndexSet>& supports)
{
    insert(supports, 0, supports.size());
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::insert(const std::vector<IndexSet>& supports, Index start, Index end)
{
    for (int i = start; i < end; ++i) { insert(supports[i], i); }
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::insert(const IndexSet& support, Index index)
{
    supp_size = support.get_size();
    root = root->insert(support, index);
}

template <class IndexSet>
inline
bool
MultiTree<IndexSet>::dominated(const IndexSet& support, Index index1, Index index2) const
{
    return root->dominated(support, index1, index2);
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::dump() const
{
    *out << "Multi Tree Dump:\n";
    root->dump(0);
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    root->find_singleton_diff(inds, s);
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    root->find(inds, zeros, s, count);
}

template <class IndexSet>
inline
void
MultiTree<IndexSet>::find(std::vector<std::pair<Index,Index> >& inds, const MultiTree& tree, Size count) const
{
    IndexSet s(supp_size,0);
    root->find(inds, s, tree.root, count);
}

////////////////
// TreeBranch //
////////////////

template <class IndexSet>
MultiTree<IndexSet>::TreeBranch::TreeBranch(Index index, const IndexSet& s1, TreeNode* n1, const IndexSet& s0, TreeNode* n0)
{
    indices.push_back(index);
    supps.push_back(s1);
    nodes.push_back(n1);
    supps.push_back(s0);
    nodes.push_back(n0);
}

template <class IndexSet>
MultiTree<IndexSet>::TreeBranch::~TreeBranch()
{
    for (size_t i = 0; i < nodes.size(); ++i)  { delete nodes[i]; }
    indices.clear();
    nodes.clear();
    supps.clear();
}

template <class IndexSet>
typename MultiTree<IndexSet>::TreeNode*
MultiTree<IndexSet>::TreeBranch::insert(const IndexSet& s, Index index)
{
    size_t i = 0;
    while (i < indices.size() && !s[indices[i]]) { ++i; }
    nodes[i] = nodes[i]->insert(s, index);
    supps[i].set_intersection(s);

#if 1
    // Collapse branches.
    if (i == indices.size()) {
        TreeBranch* branch = dynamic_cast<TreeBranch*>(nodes[i]);
        if (branch) {
            nodes.pop_back();
            supps.pop_back();
            nodes.insert(nodes.end(), branch->nodes.begin(), branch->nodes.end());
            branch->nodes.clear();
            supps.insert(supps.end(), branch->supps.begin(), branch->supps.end());
            branch->supps.clear();
            indices.insert(indices.end(), branch->indices.begin(), branch->indices.end());
            branch->indices.clear();
        }
    }
#endif
    return this;
}

template <class IndexSet>
bool
MultiTree<IndexSet>::TreeBranch::dominated(const IndexSet& s, Index index1, Index index2) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_subset(s) && nodes[i]->dominated(s, index1, index2)) { return true; }
    }
    return false;
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBranch::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].singleton_diff(s)) { nodes[i]->find_singleton_diff(inds, s); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBranch::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_disjoint(zeros) && supps[i].count_union(s) <= count) {
            nodes[i]->find(inds, zeros, s, count);
        }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBranch::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, const TreeNode* tree, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].count_union(s) <= count) { tree->find(inds, supps[i], nodes[i], count); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBranch::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, Index index, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].count_union(s) <= count) { nodes[i]->find(inds, s, index, count); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBranch::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Limb:\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (i < indices.size()) { *out << indent << "-" << indices[i] << ":"; }
        else { *out << indent << "--" << ":"; }
        for (Index j = 0; j < supps[i].get_size(); ++j) {
            if (supps[i][j]) { *out << " " << j; }
        }
        *out << "\n";
        nodes[i]->dump(level+1);
    }
}

//////////////////////
// TreeBucket //
//////////////////////

template <class IndexSet>
Size MultiTree<IndexSet>::TreeBucket::MAX_BUCKET_SIZE = 50;

template <class IndexSet>
MultiTree<IndexSet>::TreeBucket::TreeBucket()
{
}

template <class IndexSet>
MultiTree<IndexSet>::TreeBucket::~TreeBucket()
{
}

template <class IndexSet>
typename MultiTree<IndexSet>::TreeNode*
MultiTree<IndexSet>::TreeBucket::insert(const IndexSet& s, Index index)
{
    supps.push_back(s);
    indices.push_back(index);

    if ((Size) supps.size() >= MAX_BUCKET_SIZE) {
        int count = -1;
        int goal = supps.size()/2;
        Index new_index = -1;
        for (Index i = 0; i < s.get_size(); ++i) {
            int tmp_count = 0;
            for (Index j = 0; j < (Index) supps.size(); ++j) { if (supps[j][i]) { ++tmp_count; } }
            if (abs(goal-tmp_count) < abs(goal-count)) { new_index = i; count = tmp_count; }
        }

        if (count == 0 || count == (Index) supps.size()) { return this; }

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
        TreeBranch* new_branch = new TreeBranch(new_index, is1, one_bucket, is0, this);
        return new_branch;
    }
    return this;
}

template <class IndexSet>
bool
MultiTree<IndexSet>::TreeBucket::dominated(const IndexSet& s, Index index1, Index index2) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_subset(s) && indices[i] != index1 && indices[i] != index2) { return true; }
    }
    return false;
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBucket::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].singleton_diff(s)) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBucket::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_disjoint(zeros) && supps[i].count_union(s) <= count) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBucket::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, const TreeNode* tree, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].count_union(s) <= count) { tree->find(inds, supps[i], indices[i], count); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBucket::find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& s, Index index, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].count_union(s) <= count) { inds.push_back(std::pair<Index,Index>(indices[i],index)); }
    }
}

template <class IndexSet>
void
MultiTree<IndexSet>::TreeBucket::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        *out << indent << "[" << indices[i] << "] " << supps[i] << "\n";
    }
}
