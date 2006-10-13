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

#include "OnesTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

template <class IndexSet>
const int OnesTreeNode<IndexSet>::INDENT = 4;

template <class IndexSet>
OnesTreeNode<IndexSet>::OnesTreeNode()
{
}

template <class IndexSet>
OnesTreeNode<IndexSet>::~OnesTreeNode()
{
}

template <class IndexSet>
OnesTreeBranch<IndexSet>::OnesTreeBranch()
        : OnesTreeNode<IndexSet>()
{
}

template <class IndexSet>
OnesTreeBranch<IndexSet>::~OnesTreeBranch()
{
    for (unsigned int i = 0; i < nodes.size(); ++i)  { delete nodes[i].second; }
}

template <class IndexSet>
void
OnesTreeBranch<IndexSet>::insert(const IndexSet& support, int start, int remaining, int index)
{
    assert(remaining != 0);
    // There should always be another unless remaining == 0.
    int next_one = start;
    while (!support[next_one]) { ++next_one; }

    int i = 0;
    while (i < (int) nodes.size() && next_one != nodes[i].first) { ++i; }
    if (i < (int) nodes.size())
    {
        nodes[i].second->insert(support, next_one+1, remaining-1, index);
    }
    else
    {
        OnesTreeNode<IndexSet>* new_node;
        if (remaining > 1) { new_node = new OnesTreeBranch<IndexSet>(); }
        else { new_node = new OnesTreeLeaf<IndexSet>(); }
        nodes.push_back(std::pair<int,OnesTreeNode<IndexSet>*>(next_one,new_node));
        new_node->insert(support, next_one+1, remaining-1, index);
    }
}

template <class IndexSet>
bool
OnesTreeBranch<IndexSet>::dominated(const IndexSet& b, int index1, int index2)
{
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if (b[nodes[i].first] && nodes[i].second->dominated(b, index1, index2)) { return true; }
    }
    return false;
}

template <class IndexSet>
int
OnesTreeBranch<IndexSet>::index_dominated(const IndexSet& b, int index1, int index2)
{
    int index;
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if (b[nodes[i].first])
        {
            index = nodes[i].second->index_dominated(b, index1, index2);
            if (index >= 0) { return index; }
        }
    }
    return -1;
}

template <class IndexSet>
int
OnesTreeBranch<IndexSet>::find(const IndexSet& support, int start, int remaining)
{
    assert(remaining != 0);
    // There should always be another unless remaining == 0.
    int next_one = start;
    while (!support[next_one]) { ++next_one; }

    int i = 0;
    while (i < (int) nodes.size() && next_one != nodes[i].first) { ++i; }
    if (i < (int) nodes.size())
    {
        return nodes[i].second->find(support, next_one+1, remaining-1);
    }
    return -1;
}

template <class IndexSet>
void
OnesTreeBranch<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp, int diff)
{
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if (supp[nodes[i].first])
        {
            if (diff > 0) { nodes[i].second->find_diff(indices, supp, diff-1); }
        }
        else
        {
            nodes[i].second->find_diff(indices, supp, diff);
        }
    }
}

template <class IndexSet>
void
OnesTreeBranch<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2)
{
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
#if 0
        int temp_diff1 = diff1;
        int temp_diff2 = diff2;
        if (supp1[nodes[i].first]) { --temp_diff1; }
        if (temp_diff1 < 0) { continue; }
        if (supp2[nodes[i].first]) { --temp_diff2; }
        if (temp_diff2 < 0) { continue; }
        nodes[i].second->find_diff(indices, supp1, temp_diff1, supp2, temp_diff2);
#endif
#if 1
        if (supp1[nodes[i].first])
        {
            if (diff1 > 0)
            {
                if (supp2[nodes[i].first])
                {
                    if (diff2 > 0)
                    {
                        nodes[i].second->find_diff(indices, supp1, diff1-1, supp2, diff2-1);
                    }
                }
                else
                {
                    nodes[i].second->find_diff(indices, supp1, diff1-1, supp2, diff2);
                }
            }
        }
        else
        {
            if (supp2[nodes[i].first])
            {
                if (diff2 > 0)
                {
                    nodes[i].second->find_diff(indices, supp1, diff1, supp2, diff2-1);
                }
            }
            else
            {
                nodes[i].second->find_diff(indices, supp1, diff1, supp2, diff2);
            }
        }
#endif
    }
}

template <class IndexSet>
void
OnesTreeBranch<IndexSet>::dump(int level)
{
    std::string indent(level*OnesTreeNode<IndexSet>::INDENT, ' ');
    //*out << indent << level << "\n";
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        *out << indent << nodes[i].first << "\n";
        nodes[i].second->dump(level+1);
    }
}

template <class IndexSet>
OnesTreeLeaf<IndexSet>::OnesTreeLeaf()
        : OnesTreeNode<IndexSet>()
{
    index = -1;
}

template <class IndexSet>
OnesTreeLeaf<IndexSet>::~OnesTreeLeaf()
{
}

template <class IndexSet>
void
OnesTreeLeaf<IndexSet>::insert(const IndexSet& support, int start, int remaining, int _index)
{
    // If index != -1, then this implies that the leaf already exists.
    assert(index == -1);
    assert(remaining == 0);
    index = _index;
    //*out << "Start = " << start << " Rem = " << remaining;
    //*out << " Index = " << index << "\n"; 
}

template <class IndexSet>
bool
OnesTreeLeaf<IndexSet>::dominated(const IndexSet& b, int index1, int index2)
{
    if (index != index1 && index != index2) { return true; }
    return false;
}

template <class IndexSet>
int
OnesTreeLeaf<IndexSet>::index_dominated(const IndexSet& b, int index1, int index2)
{
    if (index != index1 && index != index2) { return index; }
    return -1;
}

template <class IndexSet>
int
OnesTreeLeaf<IndexSet>::find(const IndexSet& support, int start, int remaining)
{
    if (remaining == 0) { return index; }
    return -1;
}

template <class IndexSet>
void
OnesTreeLeaf<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp, int diff)
{
    indices.push_back(index);
}

template <class IndexSet>
void
OnesTreeLeaf<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2)
{
    indices.push_back(index);
}

template <class IndexSet>
void
OnesTreeLeaf<IndexSet>::dump(int level)
{
    std::string indent(level*OnesTreeNode<IndexSet>::INDENT, ' ');
    *out << indent << index << "\n";
}

template <class IndexSet>
OnesTree<IndexSet>::OnesTree(const std::vector<IndexSet>& supports, int num)
{
    root = new OnesTreeBranch<IndexSet>();
    insert(supports,num);
}

template <class IndexSet>
void
OnesTree<IndexSet>::insert(const std::vector<IndexSet>& supports, int num)
{
    assert(!supports.empty());
    for (int i = 0; i < num; ++i)
    {
        //*out << "Inserting Support:\n" << supports[i] << "\n";
        root->insert(supports[i], 0, supports[i].count(), i);
    }
}

template <class IndexSet>
OnesTree<IndexSet>::OnesTree()
{
    root = new OnesTreeBranch<IndexSet>();
}

template <class IndexSet>
OnesTree<IndexSet>::~OnesTree()
{
    delete root;
}

template <class IndexSet>
bool
OnesTree<IndexSet>::dominated(const IndexSet& b, int index1, int index2)
{
    return root->dominated(b, index1, index2);
}

template <class IndexSet>
int
OnesTree<IndexSet>::index_dominated(const IndexSet& b, int index1, int index2)
{
    return root->index_dominated(b, index1, index2);
}

template <class IndexSet>
void
OnesTree<IndexSet>::insert(const IndexSet& support, int index)
{
    root->insert(support, 0, support.count(), index);
}

template <class IndexSet>
void
OnesTree<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp, int diff)
{
    root->find_diff(indices, supp, diff);
}

template <class IndexSet>
void
OnesTree<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2)
{
    root->find_diff(indices, supp1, diff1, supp2, diff2);
}

template <class IndexSet>
void
OnesTree<IndexSet>::dump()
{
    *out << "Tree Dump:\n";
    root->dump(0);
}


