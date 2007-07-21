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

#include "SupportTree.h"

#include <iostream>
#include <string>

using namespace _4ti2_;

template <class IndexSet>
const int SupportTree<IndexSet>::INDENT = 4;

template <class IndexSet>
SupportTree<IndexSet>::SupportTreeNode::SupportTreeNode()
{
    index = -1;
}

template <class IndexSet>
SupportTree<IndexSet>::SupportTreeNode::~SupportTreeNode()
{
    for (unsigned int i = 0; i < nodes.size(); ++i)  { delete nodes[i].second; }
}

template <class IndexSet>
void
SupportTree<IndexSet>::insert(
                SupportTreeNode& node, const IndexSet& support,
                int start, int remaining, int index)
{
    assert(index >= 0);
    if (remaining > 0)
    {
        // There should always be another unless remaining == 0.
        int next_one = start;
        while (!support[next_one]) { ++next_one; }
    
        int i = 0;
        while (i < (int) node.nodes.size() && next_one != node.nodes[i].first) { ++i; }
        if (i < (int) node.nodes.size())
        {
            insert(*node.nodes[i].second, support, next_one+1, remaining-1, index);
        }
        else
        {
            SupportTreeNode* new_node = new SupportTreeNode;
            node.nodes.push_back(std::pair<int,SupportTreeNode*>(next_one,new_node));
            insert(*new_node, support, next_one+1, remaining-1, index);
        }
    }
    else
    {
        assert(node.index == -1);
        node.index = index;
    }
}

template <class IndexSet>
bool
SupportTree<IndexSet>::dominated(SupportTreeNode& node, const IndexSet& b, int index1, int index2)
{
    if (node.index < 0)
    {
        for (unsigned int i = 0; i < node.nodes.size(); ++i)
        {
            if (b[node.nodes[i].first] && dominated(*node.nodes[i].second, b, index1, index2))
            { return true; }
        }
        return false;
    }
    else
    {
        if (node.index != index1 && node.index != index2) { return true; }
        return false;
    }
}

template <class IndexSet>
int
SupportTree<IndexSet>::index_dominated(SupportTreeNode& node, const IndexSet& b, int index1, int index2)
{
    if (node.index < 0)
    {
        int index;
        for (unsigned int i = 0; i < node.nodes.size(); ++i)
        {
            if (b[node.nodes[i].first])
            {
                index = index_dominated(*node.nodes[i].second, b, index1, index2);
                if (node.index >= 0) { return index; }
            }
        }
        return -1;
    }
    else
    {
        if (node.index != index1 && node.index != index2) { return node.index; }
        return -1;
    }
}

template <class IndexSet>
int
SupportTree<IndexSet>::find(SupportTreeNode& node, const IndexSet& support, int start, int remaining)
{
    if (remaining > 0)
    {
        // There should always be another unless remaining == 0.
        int next_one = start;
        while (!support[next_one]) { ++next_one; }
    
        int i = 0;
        while (i < (int) node.nodes.size() && next_one != node.nodes[i].first) { ++i; }
        if (i < (int) node.nodes.size())
        {
            return find(*node.nodes[i].second, support, next_one+1, remaining-1);
        }
        return -1;
    }
    else
    {
        return node.index;
    }
}

template <class IndexSet>
void
SupportTree<IndexSet>::find_diff(
                SupportTreeNode& node, std::vector<int>& indices,
                const IndexSet& supp, int diff)
{
    if (node.index < 0)
    {
        for (unsigned int i = 0; i < node.nodes.size(); ++i)
        {
            if (supp[node.nodes[i].first])
            {
                if (diff > 0) { find_diff(*node.nodes[i].second, indices, supp, diff-1); }
            }
            else
            {
                find_diff(*node.nodes[i].second, indices, supp, diff);
            }
        }
    }
    else
    {
        indices.push_back(node.index);
    }
}

template <class IndexSet>
void
SupportTree<IndexSet>::find_diff(
                SupportTreeNode& node, std::vector<int>& indices,
                const IndexSet& supp1, int diff1,
                const IndexSet& supp2, int diff2)
{
    if (node.index < 0)
    {
        for (unsigned int i = 0; i < node.nodes.size(); ++i)
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
            if (supp1[node.nodes[i].first])
            {
                if (diff1 > 0)
                {
                    if (supp2[node.nodes[i].first])
                    {
                        if (diff2 > 0)
                        {
                            find_diff(*node.nodes[i].second, indices, supp1, diff1-1, supp2, diff2-1);
                        }
                    }
                    else
                    {
                        find_diff(*node.nodes[i].second, indices, supp1, diff1-1, supp2, diff2);
                    }
                }
            }
            else
            {
                if (supp2[node.nodes[i].first])
                {
                    if (diff2 > 0)
                    {
                        find_diff(*node.nodes[i].second, indices, supp1, diff1, supp2, diff2-1);
                    }
                }
                else
                {
                    find_diff(*node.nodes[i].second, indices, supp1, diff1, supp2, diff2);
                }
            }
#endif
        }
    }
    else
    {
        indices.push_back(node.index);
    }
}

template <class IndexSet>
void
SupportTree<IndexSet>::dump(SupportTreeNode& node, int level)
{
    std::string indent(level*INDENT, ' ');
    if (node.index < 0)
    {
        assert(node.nodes.size() > 0);
        for (unsigned int i = 0; i < node.nodes.size(); ++i)
        {
            *out << indent << node.nodes[i].first << "\n";
            dump(*node.nodes[i].second, level+1);
        }
    }
    else
    {
        *out << indent << "index = " << node.index << "\n";
    }
}

template <class IndexSet>
SupportTree<IndexSet>::SupportTree(const std::vector<IndexSet>& supports, int num)
{
    insert(supports, num);
}


template <class IndexSet>
SupportTree<IndexSet>::SupportTree()
{
}

template <class IndexSet>
SupportTree<IndexSet>::~SupportTree()
{
}
