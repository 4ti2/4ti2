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

#ifndef _4ti2_qsolve__BinaryTree_
#define _4ti2_qsolve__BinaryTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class BinaryTree
{
public:
    BinaryTree(const std::vector<IndexSet>& supports, int num);
    ~BinaryTree();

    void insert(const std::vector<IndexSet>& supports, int num);
    bool dominated(const IndexSet& b, int index1, int index2);
    void dump();
    void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff);

private:
    class BinaryTreeNode
    {
    public:
        BinaryTreeNode() {}
        virtual ~BinaryTreeNode() {}
        virtual BinaryTreeNode* insert(const IndexSet& support, int index) = 0;
        virtual bool dominated(const IndexSet& b, int index1, int index2) = 0;
        virtual void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff) = 0;
        virtual void dump(int level) = 0;
    };

    class BinaryTreeLeaf : public BinaryTreeNode {
    public:
        BinaryTreeLeaf(Index _index, const IndexSet& _is) : i(_index), is(_is) {}
        ~BinaryTreeLeaf() {}
        BinaryTreeNode* insert(const IndexSet& support, int index);
        bool dominated(const IndexSet& b, int index1, int index2) {
            if (IndexSet::set_subset(is, b) && i != index1 && i != index2) { return true; }
            return false;
        }
        void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff) {
            IndexSet tmp(is);
            tmp.set_intersection(supp);
            if (tmp.count() <= orig_diff) { indices.push_back(i); }
        }
        void dump(int level) {
            std::string indent(level*INDENT, ' ');
            *out << indent << "[" << i << "] " << is << "\n";
        }
    private:
        Index i;
        IndexSet is;
    };

    class BinaryTreeBranch : public BinaryTreeNode {
    public:
        BinaryTreeBranch(Index index, BinaryTreeNode* _zero, BinaryTreeNode* _one);
        ~BinaryTreeBranch();
        BinaryTreeNode* insert(const IndexSet& support, int index);
        bool dominated(const IndexSet& b, int index1, int index2);
        void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff);
        void dump(int level) {
            std::string indent(level*INDENT, ' ');
            *out << indent << i << " = 0:\n";
            zero->dump(level+1);
            *out << indent << i << " = 1:\n";
            one->dump(level+1);
        }
    private:
        Index i;
        BinaryTreeNode* zero;
        BinaryTreeNode* one;
    };

    void insert(const IndexSet& support, int index);

    static const int INDENT;

    BinaryTreeNode* root;
};

template <class IndexSet>
BinaryTree<IndexSet>::BinaryTree(const std::vector<IndexSet>& supports, int num)
        : root(0)
{
    insert(supports, num);
}

template <class IndexSet>
BinaryTree<IndexSet>::~BinaryTree()
{
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::insert(const std::vector<IndexSet>& supports, int num)
{
    for (int i = 0; i < num; ++i) {
        insert(supports[i], i);
    }
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::insert(const IndexSet& support, int index)
{
    if (root == 0) { root = new BinaryTreeLeaf(index, support); }
    else { root = root->insert(support, index); }
}

template <class IndexSet>
inline
bool
BinaryTree<IndexSet>::dominated(const IndexSet& support, int index1, int index2)
{
    return root->dominated(support, index1, index2);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::dump()
{
    *out << "Binary Tree Dump:\n";
    root->dump(0);
}

template <class IndexSet>
inline
void
BinaryTree<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp, int diff)
{
    root->find_diff(indices, supp, diff, diff);
}

} // namespace _4ti2_

// Include template defitinions.
#include "qsolve/BinaryTree.hpp"

#endif
