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

#ifndef _4ti2_groebner__SupportTree_
#define _4ti2_groebner__SupportTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class SupportTree
{
public:
    SupportTree(const std::vector<IndexSet>& supports, int num);
    SupportTree();
    ~SupportTree();

    void insert(const std::vector<IndexSet>& supports, int num);
    void insert(const IndexSet& support, int index);
    bool dominated(const IndexSet& b);
    bool dominated(const IndexSet& b, int index1, int index2);
    int index_dominated(const IndexSet& b, int index1, int index2);
    int find(const IndexSet& b);
    void dump();

    void find_diff(std::vector<int>& indices,
                    const IndexSet& supp, int diff);
    void find_diff(std::vector<int>& indices,
                    const IndexSet& supp1, int diff1,
                    const IndexSet& supp2, int diff2);

private:
    class SupportTreeNode
    {
    public:
        SupportTreeNode();
        ~SupportTreeNode();
        std::vector<std::pair<int,SupportTreeNode*> > nodes;
        int index;
    };

    void insert(SupportTreeNode& node, const IndexSet& support, int start, int remaining, int index);
    void dump(SupportTreeNode& node, int level);
    bool dominated(SupportTreeNode& node, const IndexSet& b);
    bool dominated(SupportTreeNode& node, const IndexSet& b, int index1, int index2);
    int index_dominated(SupportTreeNode& node, const IndexSet& b, int index1, int index2);
    int find(SupportTreeNode& node, const IndexSet& b, int start, int remaining);

    void find_diff(SupportTreeNode& node, std::vector<int>& indices,
                    const IndexSet& supp, int diff);
    void find_diff(SupportTreeNode& node, std::vector<int>& indices,
                    const IndexSet& supp1, int diff1,
                    const IndexSet& supp2, int diff2);

    static const int INDENT;

    SupportTreeNode root;
};

template <class IndexSet>
inline
void
SupportTree<IndexSet>::insert(const std::vector<IndexSet>& supports, int num)
{
    for (int i = 0; i < num; ++i)
    {
        insert(root, supports[i], 0, supports[i].count(), i);
    }
}

template <class IndexSet>
inline
void
SupportTree<IndexSet>::insert(const IndexSet& support, int index)
{
    insert(root, support, 0, support.count(), index);
}

template <class IndexSet>
inline
bool
SupportTree<IndexSet>::dominated(const IndexSet& support)
{
    return dominated(root, support);
}

template <class IndexSet>
inline
bool
SupportTree<IndexSet>::dominated(const IndexSet& support, int index1, int index2)
{
    return dominated(root, support, index1, index2);
}

template <class IndexSet>
inline
int
SupportTree<IndexSet>::index_dominated(const IndexSet& support, int index1, int index2)
{
    return index_dominated(root, support, index1, index2);
}

template <class IndexSet>
inline
int
SupportTree<IndexSet>::find(const IndexSet& support)
{
    return find(root, support, 0, support.count());
}

template <class IndexSet>
inline
void
SupportTree<IndexSet>::dump()
{
    *out << "Support Tree Dump:\n";
    dump(root, 0);
}

template <class IndexSet>
inline
void
SupportTree<IndexSet>::find_diff(
                std::vector<int>& indices,
                const IndexSet& supp, int diff)
{
    find_diff(root, indices, supp, diff);
}

template <class IndexSet>
inline
void
SupportTree<IndexSet>::find_diff(
                std::vector<int>& indices,
                const IndexSet& supp1, int diff1,
                const IndexSet& supp2, int diff2)
{
    find_diff(root, indices, supp1, diff1, supp2, diff2);
}

} // namespace _4ti2_

// Include template defitinions.
#include "groebner/SupportTree.tpp"

#endif
