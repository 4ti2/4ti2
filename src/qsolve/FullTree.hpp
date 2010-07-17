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

#include "qsolve/FullTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

template <class IndexSet>
const int FullTree<IndexSet>::INDENT = 2;
template <class IndexSet>
long long int FullTree<IndexSet>::num_nodes = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::num_buckets = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::num_limbs = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::num_buds = 0;

//////////////
// FullTree //
//////////////

template <class IndexSet>
FullTree<IndexSet>::FullTree()
{
}

template <class IndexSet>
FullTree<IndexSet>::~FullTree()
{
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::insert(const std::vector<IndexSet>& supports)
{
    for (size_t i = 0; i < supports.size(); ++i) {
        insert(supports[i], i);
    }
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::insert(const std::vector<IndexSet>& supports, Index start, Index end)
{
    for (Index i = start; i < end; ++i) {
        insert(supports[i], i);
    }
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::insert(const IndexSet& support, int index)
{
    IndexSet empty(support.get_size(),0);
    root.insert(support, empty, index);
}

template <class IndexSet>
inline
bool
FullTree<IndexSet>::dominated(const IndexSet& support, int index1, int index2) const
{
    return root.dominated(support, index1, index2);
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::dump() const
{
    *out << "Full Tree Dump:\n";
    root.dump(0);
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::clear()
{
    root.clear();
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::find_singleton_diff(std::vector<int>& inds, const IndexSet& supp) const
{
    root.find_singleton_diff(inds, supp);
}

template <class IndexSet>
inline
void
FullTree<IndexSet>::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const
{
    root.find(inds, zeros, supp, count);
}

template <class IndexSet>
void
FullTree<IndexSet>::print_statistics()
{
    num_nodes = num_limbs = num_buckets = num_buds = 0;
    root.compute_statistics();
    *out << "\nSTATISTICS:\n";
    *out << "Num Nodes = " << num_nodes << "\n";
    *out << "Num Limbs = " << num_limbs << "\n";
    *out << "Num Buckets = " << num_buckets << "\n";
    *out << "Num Buds = " << num_buds << "\n";
}

/////////////
// TreeBud //
/////////////

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeBud::insert(const IndexSet& supp, const IndexSet& parent, Index index)
{
    assert(index >= 0);

    TreeBucket* bucket = new TreeBucket();
    bucket->insert(parent, parent, i);
    bucket->insert(supp, supp, index);

    delete this;
    return bucket;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeBud::dominated(const IndexSet& b, Index index1, Index index2) const
{
    if (i != index1 && i != index2) { return true; }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const
{
    inds.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::find(std::vector<Index>& inds, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    inds.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::dump(Index level) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "B[" << i << "]\n";
}

//////////////
// TreeLimb //
//////////////

template <class IndexSet>
FullTree<IndexSet>::TreeLimb::TreeLimb()
{
}

template <class IndexSet>
FullTree<IndexSet>::TreeLimb::~TreeLimb()
{
    clear();
}

template <class IndexSet>
inline void 
FullTree<IndexSet>::TreeLimb::clear()
{
    for (size_t i = 0; i < nodes.size(); ++i)  { delete nodes[i]; }
    indices.clear();
    nodes.clear();
    supps.clear();
}

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeLimb::insert(const IndexSet& supp, const IndexSet& parent, Index index)
{
    size_t i = 0;
    while (i < indices.size() && !supp[indices[i]]) { ++i; }
    if (i == indices.size()) {
        IndexSet temp_supp(supp);
        temp_supp.set_difference(parent);
        Index index = 0;
        while (!temp_supp[index]) { ++index; }
        TreeNode* node = new TreeBud(index);
        nodes.push_back(node);
        indices.push_back(index);
        supps.push_back(supp);
    }
    else {
        nodes[i] = nodes[i]->insert(supp, supps[i], index);
        supps[i].set_intersection(supp);
    }
    return this;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::insert(const IndexSet& supp, Index index, TreeNode* node)
{
    supps.push_back(supp);
    nodes.push_back(node);
    indices.push_back(index);
    return;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeLimb::dominated(const IndexSet& s, int index1, int index2) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_subset(s) && nodes[i]->dominated(s, index1, index2)) {
            return true;
        }
    }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].singleton_diff(supp)) { nodes[i]->find_singleton_diff(inds, supp); }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::find(std::vector<int>& inds, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_disjoint(zero) && supps[i].count_union(supp) <= count) {
            nodes[i]->find(inds, zero, supp, count);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::dump(Index level) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Limb:\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        *out << indent << "-" << indices[i] << ":";
        for (Index j = 0; j < supps[i].get_size(); ++j) {
            if (supps[i][j]) { *out << " " << j; }
        }
        *out << "\n";
        nodes[i]->dump(level+1);
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::compute_statistics()
{
    ++num_nodes;
    ++num_limbs;
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i]->compute_statistics();
    }
}

////////////////
// TreeBucket //
////////////////

template <class IndexSet>
const Size FullTree<IndexSet>::TreeBucket::MAX_BUCKET_SIZE = 50;

template <class IndexSet>
FullTree<IndexSet>::TreeBucket::TreeBucket()
{
}

template <class IndexSet>
FullTree<IndexSet>::TreeBucket::~TreeBucket()
{
}

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeBucket::insert(const IndexSet& s, const IndexSet& p, Index index)
{
    supps.push_back(s);
    indices.push_back(index);

    if ((Size) supps.size() >= MAX_BUCKET_SIZE) {
        TreeLimb* limb = new TreeLimb();

        while (!supps.empty()) {
            int count = -1;
            int goal = supps.size()/2;
            Index new_index = -1;
            for (Index i = 0; i < s.get_size(); ++i) {
                int tmp_count = 0;
                for (Index j = 0; j < (Index) supps.size(); ++j) { if (supps[j][i]) { ++tmp_count; } }
                if (abs(goal-tmp_count) < abs(goal-count)) { new_index = i; count = tmp_count; }
            }

            if (count == 1) {
                for (Index i = supps.size()-1; i >= 0; --i) {
                    if (supps[i][new_index]) {
                        TreeBud* bud = new TreeBud(indices[i]);
                        limb->insert(supps[i], new_index, bud);
                        supps.erase(supps.begin()+i);
                        indices.erase(indices.begin()+i);
                        break;
                    }
                }
            }
            else {
                IndexSet is(s.get_size(),1);
                TreeBucket* bucket = new TreeBucket();
                for (Index i = supps.size()-1; i >= 0; --i) {
                    if (supps[i][new_index]) {
                        is.set_intersection(supps[i]);
                        bucket->insert(supps[i], supps[i], indices[i]);
                        supps.erase(supps.begin()+i);
                        indices.erase(indices.begin()+i);
                    }
                }
                limb->insert(is, new_index, bucket);
            }
        }
        delete this;
        return limb;
    }
    return this;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeBucket::dominated(const IndexSet& s, Index index1, Index index2) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_subset(s) && indices[i] != index1 && indices[i] != index2) { return true; }
    }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBucket::find_singleton_diff(std::vector<Index>& inds, const IndexSet& s) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].singleton_diff(s)) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBucket::find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& s, Size count) const
{
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        if (supps[i].set_disjoint(zeros) && supps[i].count_union(s) <= count) { inds.push_back(indices[i]); }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBucket::dump(int level) const
{
    std::string indent(level*INDENT, ' ');
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        *out << indent << "[" << indices[i] << "] " << supps[i] << "\n";
    }
}
