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
long long int FullTree<IndexSet>::TreeNode::num_nodes = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_branches = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_limbs = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_arms = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_leaves = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_buds = 0;
template <class IndexSet>
long long int FullTree<IndexSet>::TreeNode::num_singletons = 0;

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
    //if (index == 1145) { 
    //    *out << "\n";   
    //    *out << "Adding[" << index << "]: ";
    //    for (Index i = 0; i < support.get_size(); ++i) {
    //        if (support[i]) { *out << i << " "; }
    //    }
    //    *out << "\n";   
    //    dump();
    //}
    root.insert(support, 0, index);
    //if (index == 1145) { dump(); }
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
    root.dump(0, 0);
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
FullTree<IndexSet>::find_diff(std::vector<int>& inds, const IndexSet& supp, int diff) const
{
    root.find_diff(inds, supp, diff, diff);
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
inline
void
FullTree<IndexSet>::find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const
{
    root.find(inds, inds0, zeros, supp, count);
}

template <class IndexSet>
void
FullTree<IndexSet>::print_statistics()
{
    root.reset_statistics();
    root.compute_statistics();
    root.print_statistics();
}

//////////////
// TreeNode //
//////////////

template <class IndexSet>
void
FullTree<IndexSet>::TreeNode::reset_statistics()
{
    num_nodes = num_branches = num_limbs = num_arms = num_leaves = num_buds = num_singletons = 0;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeNode::print_statistics()
{
    *out << "\nSTATISTICS:\n";
    *out << "Num Nodes = " << num_nodes << "\n";
    *out << "Num Branches = " << num_branches << "\n";
    *out << "Num Limbs = " << num_limbs << "\n";
    *out << "Num Arms = " << num_arms << "\n";
    *out << "Num Leaves = " << num_leaves << "\n";
    *out << "Num Buds = " << num_buds << "\n";
    *out << "Num Singletons = " << num_singletons << "\n";
}

////////////////
// TreeBranch //
////////////////

#if 0
template <class IndexSet>
inline
FullTree<IndexSet>::TreeBranch::TreeBranch()
{
}

template <class IndexSet>
inline
FullTree<IndexSet>::TreeBranch::TreeBranch(Index next, TreeNode* node)
{
    nodes.push_back(std::pair<Index,TreeNode*>(next,node));
}

template <class IndexSet>
inline
FullTree<IndexSet>::TreeBranch::~TreeBranch()
{
    clear();
}

template <class IndexSet>
inline void 
FullTree<IndexSet>::TreeBranch::clear()
{
    for (size_t i = 0; i < nodes.size(); ++i)  { delete nodes[i].second; }
    nodes.clear();
}

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeBranch::insert(const IndexSet& support, int start, int index)
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
        TreeNode* new_node = new TreeBud(index, support);
        nodes.insert(nodes.begin()+i, std::pair<Index,TreeNode*>(next_one,new_node));
    }
    return this;

#if 0
    size_t i = 0;
    while (i < nodes.size() && next_one != nodes[i].first) { ++i; }
    if (i < nodes.size()) {
        nodes[i].second = nodes[i].second->insert(support, next_one+1, index);
    }
    else {
        //TreeNode* new_node = new TreeLeaf(index, support);
        TreeNode* new_node = new TreeBud(index, support);
        nodes.push_back(std::pair<Index,TreeNode*>(next_one,new_node));
    }
    return this;
#endif
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeBranch::dominated(const IndexSet& b, int index1, int index2) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (b[nodes[i].first] && nodes[i].second->dominated(b, index1, index2)) {
            return true;
        }
    }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBranch::find_diff(
            std::vector<int>& inds, const IndexSet& supp, int diff, int orig_diff) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supp[nodes[i].first]) {
            if (diff > 0) { nodes[i].second->find_diff(inds, supp, diff-1, orig_diff); }
        }
        else {
            nodes[i].second->find_diff(inds, supp, diff, orig_diff);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBranch::find_singleton_diff(std::vector<int>& inds, const IndexSet& supp) const
{
    *out << "ERROR: function `find_singleton_diff' not support for TreeBranch." << std::endl;
    exit(1);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBranch::dump(Index level, Index index) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Branch:\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        *out << indent << "-" << nodes[i].first << "\n";
        nodes[i].second->dump(level+1, nodes[i].first+1);
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBranch::compute_statistics()
{
    ++TreeNode::num_nodes;
    ++TreeNode::num_branches;
    if (nodes.size() == 1) { ++TreeNode::num_singletons; }
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].second->compute_statistics();
    }
}
#endif

//////////////
// TreeLeaf //
//////////////

#if 0
template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeLeaf::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);

    // There should always be another.
    int next_one = start;
    while (!is[next_one]) { ++next_one; }

    TreeBranch* new_branch = new TreeBranch(next_one, this);
    new_branch->insert(support, start, index);

    return new_branch;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeLeaf::dominated(const IndexSet& b, Index index1, Index index2) const
{
    if (is.set_subset(b) && i != index1 && i != index2) { return true; }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLeaf::find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const 
{
    IndexSet tmp(is);
    tmp.set_intersection(supp);
    if (tmp.singleton()) { inds.push_back(i); }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLeaf::find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const
{
    if (is.singleton_diff(supp)) { inds.push_back(i); }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLeaf::dump(Index level, Index index) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Leaf[" << i << "]: ";
    for (Index i = index; i < is.get_size(); ++i) {
        if (is[i]) { *out << i << " "; }
    }
    *out << "\n";
}
#endif

/////////////
// TreeBud //
/////////////

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeBud::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);

    //TreeBranch* new_branch = new TreeBranch(next_one, this);
    TreeBud* new_bud = new TreeBud(index, support);
    TreeLimb* new_branch = new TreeLimb(start, support, new_bud, is, this);

    return new_branch;
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
FullTree<IndexSet>::TreeBud::find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const 
{
    IndexSet tmp(is);
    tmp.set_intersection(supp);
    if (tmp.singleton()) { inds.push_back(i); }
    //indices.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const
{
    inds.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::find(std::vector<int>& inds, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    inds.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const
{
    if (supp.singleton_diff(is)) { inds0.push_back(i);  return; }
    if (is.count_lte_diff(2, supp)) { inds0.push_back(i); return; }
    inds.push_back(i);
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeBud::dump(Index level, Index index) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "B[" << i << "]: ";
    for (Index i = index; i < is.get_size(); ++i) {
        if (is[i]) { *out << i << " "; }
    }
    *out << "\n";
}

//////////////
// TreeLimb //
//////////////

template <class IndexSet>
FullTree<IndexSet>::TreeLimb::TreeLimb(Index i,
            const IndexSet& s1, TreeNode* n1,
            const IndexSet& s2, TreeNode* n2)
    : supp_intersection(s1.get_size())
{
    supp_intersection.set_intersection(s1, s2);
    // Find the first index after i where s1 and s2 differ.
    IndexSet temp_supp(s1.get_size());
    temp_supp.set_symm_difference(s1, s2);
    while (!temp_supp[i] && i < temp_supp.get_size()) { ++i; }
    //if (i == temp_supp.get_size()) { *out << "\nHELP1\n"; }

    if (s1[i]) {
        indices.push_back(i);
        nodes.push_back(n1);
        supps.push_back(s1);
        while (!s2[i] && i < s2.get_size()) { ++i; }
        //if (i == s2.get_size()) { *out << "\nHELP2\n"; }
        indices.push_back(i);
        nodes.push_back(n2);
        supps.push_back(s2);
    }
    else {
        indices.push_back(i);
        nodes.push_back(n2);
        supps.push_back(s2);
        while (!s1[i] && i < s1.get_size()) { ++i; }
        //if (i == s1.get_size()) { *out << "\nHELP3\n"; }
        indices.push_back(i);
        nodes.push_back(n1);
        supps.push_back(s1);
    }
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
FullTree<IndexSet>::TreeLimb::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);
    assert(!indices.empty());

    // Find the first index of difference.
    Index next = start;
    IndexSet temp_supp(support.get_size());
    temp_supp.set_symm_difference(support, supp_intersection);
    while (next < indices[0] && !temp_supp[next]) { ++next; }
    while (next < indices[0] && !supp_intersection[next]) { ++next; }
    if (next < indices[0]) {
        TreeBud* new_bud = new TreeBud(index, support); 
        TreeLimb* new_limb = new TreeLimb(start, support, new_bud, supp_intersection, this);
        return new_limb;
    }

    // There should always be another.
    next = start;
    while (!support[next] || (next < indices[0] && supp_intersection[next])) { ++next; }
    size_t i = 0;
    while (i < indices.size() && next > indices[i]) { ++i; }
    if (i < indices.size() && next == indices[i]) {
        nodes[i] = nodes[i]->insert(support, next+1, index);
        supps[i].set_intersection(support);
    }
    else {
        TreeNode* new_node = new TreeBud(index, support);
        indices.insert(indices.begin()+i, next);
        nodes.insert(nodes.begin()+i, new_node);
        supps.insert(supps.begin()+i, support);
    }
    supp_intersection.set_intersection(support);
    return this;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeLimb::dominated(const IndexSet& b, int index1, int index2) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_subset(b) && nodes[i]->dominated(b, index1, index2)) {
            return true;
        }
    }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::find_diff(
            std::vector<int>& inds, const IndexSet& supp, int diff, int orig_diff) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supp[indices[i]]) {
            if (diff > 0) { nodes[i]->find_diff(inds, supp, diff-1, orig_diff); }
        }
        else {
            nodes[i]->find_diff(inds, supp, diff, orig_diff);
        }
    }
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
FullTree<IndexSet>::TreeLimb::find(std::vector<int>& inds, std::vector<int>& inds0, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_disjoint(zero) && supps[i].count_union(supp) <= count) {
            nodes[i]->find(inds, inds0, zero, supp, count);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::dump(Index level, Index index) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Limb:";
    for (Index i = index; i < indices[0]; ++i) {
        if (supp_intersection[i]) { *out << " " << i; }
    }
    //*out << indent << "Counts:" << supp_intersection.count() << " -> ";
    //for (size_t i = 0; i < nodes.size(); ++i) {
    //    *out << indices[i] << ":" << supps[i].count() << " ";
    //}
    *out << "\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        *out << indent << "-" << indices[i] << ":";
        for (Index j = 0; j < supps[i].get_size(); ++j) {
            if (supps[i][j]) { *out << " " << j; }
        }
        *out << "\n";
        nodes[i]->dump(level+1, indices[i]+1);
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeLimb::compute_statistics()
{
    ++TreeNode::num_nodes;
    ++TreeNode::num_limbs;
    if (nodes.size() == 1) { ++TreeNode::num_singletons; }
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i]->compute_statistics();
    }
}

//////////////
// TreeArm //
//////////////

template <class IndexSet>
FullTree<IndexSet>::TreeArm::TreeArm()
{
}

template <class IndexSet>
FullTree<IndexSet>::TreeArm::~TreeArm()
{
    clear();
}

template <class IndexSet>
inline void 
FullTree<IndexSet>::TreeArm::clear()
{
    for (size_t i = 0; i < nodes.size(); ++i)  { delete nodes[i]; }
    indices.clear();
    nodes.clear();
    supps.clear();
}

template <class IndexSet>
typename FullTree<IndexSet>::TreeNode*
FullTree<IndexSet>::TreeArm::insert(const IndexSet& support, int start, int index)
{
    assert(index >= 0);
    assert(!indices.empty());

    // There should always be another.
    Index next = start;
    while (!support[next]) { ++next; }
    size_t i = 0;
    while (i < indices.size() && next > indices[i]) { ++i; }
    if (i < indices.size() && next == indices[i]) {
        nodes[i] = nodes[i]->insert(support, next+1, index);
        supps[i].set_intersection(support);
    }
    else {
        TreeNode* new_node = new TreeBud(index, support);
        indices.insert(indices.begin()+i, next);
        nodes.insert(nodes.begin()+i, new_node);
        supps.insert(supps.begin()+i, support);
    }
    return this;
}

template <class IndexSet>
bool
FullTree<IndexSet>::TreeArm::dominated(const IndexSet& b, int index1, int index2) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_subset(b) && nodes[i]->dominated(b, index1, index2)) {
            return true;
        }
    }
    return false;
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::find_diff(
            std::vector<int>& inds, const IndexSet& supp, int diff, int orig_diff) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supp[indices[i]]) {
            if (diff > 0) { nodes[i]->find_diff(inds, supp, diff-1, orig_diff); }
        }
        else {
            nodes[i]->find_diff(inds, supp, diff, orig_diff);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::find_singleton_diff(std::vector<int>& inds, const IndexSet& supp) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].singleton_diff(supp)) { nodes[i]->find_singleton_diff(inds, supp); }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::find(std::vector<int>& inds, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_disjoint(zero) && supps[i].count_union(supp) <= count) {
            nodes[i]->find(inds, zero, supp, count);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::find(std::vector<int>& inds, std::vector<int>& inds0, const IndexSet& zero, const IndexSet& supp, Size count) const
{
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (supps[i].set_disjoint(zero) && supps[i].count_union(supp) <= count) {
            nodes[i]->find(inds, inds0, zero, supp, count);
        }
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::dump(Index level, Index index) const
{
    std::string indent(level*INDENT, ' ');
    *out << indent << "Arm:" << "\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        *out << indent << "-" << indices[i] << ":";
        for (Index j = 0; j < supps[i].get_size(); ++j) {
            if (supps[i][j]) { *out << " " << j; }
        }
        *out << "\n";
        nodes[i]->dump(level+1, indices[i]+1);
    }
}

template <class IndexSet>
void
FullTree<IndexSet>::TreeArm::compute_statistics()
{
    ++TreeNode::num_nodes;
    ++TreeNode::num_arms;
    if (nodes.size() == 1) { ++TreeNode::num_singletons; }
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i]->compute_statistics();
    }
}

